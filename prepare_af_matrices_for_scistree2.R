#clean the environment####
rm(list = ls(envir = parent.frame()))
gc()


#this script will read the raw allele depth matrices, clean them and export them so they can be used with the `ScisTree2.py` script.
#libraries
library(tidyverse)
library(vroom)
library(furrr)
library(ggalign)
# pick a sensible number of workers
plan(multisession, workers = parallel::detectCores() - 1)

#set parameters
NA_threshold <- 0.4

#get the data
matrices <- list.files("inputs/nucleotide_variants/allele_depth_matrices/", pattern = "ad_matrix.tsv")

#get cells_metadata (needed for spliting matrices by strain)
cells_meta <- read_delim("inputs/cell_qc/cells_meta.tsv") %>%
  column_to_rownames("rowname") %>%
  rownames_to_column("cell") %>%
  select(sample, strain, cell) %>%
  mutate(sample = tolower(sample)) %>%
  mutate(sample = gsub(" ", "_", sample)) %>%
  mutate(simple_barcode = gsub("_.*", "", cell))

#create a function to process the matrices so it can run in parallel
process_one_matrix <- function(mfile) {
  #get matrix
  original_mat <- vroom::vroom(file.path("inputs/nucleotide_variants/allele_depth_matrices/", mfile), show_col_types = FALSE, col_types = cols(.default = col_character()))
  
  # drop all-NA columns
  keep <- !apply(original_mat, 2, function(x) all(is.na(x)))
  original_mat <- original_mat[, keep, drop = FALSE]
  
  #specify to which sample this matrix comes from 
  mat_sample <- gsub("_ad_matrix.*", "", mfile)
  
  #create a key name pair for the correct barcodes
  correct_barcodes <- cells_meta %>%
    filter(sample == mat_sample) %>%
    select(simple_barcode, cell) %>%
    deframe()
  
  #keep only cells found in the cells_meta
  keep_cols <- colnames(original_mat) %in% c("chromosome", "position", names(correct_barcodes))
  original_mat <- original_mat[, keep_cols, drop = FALSE]
  
  #correct colnames
  colnames(original_mat)[3:ncol(original_mat)] <- correct_barcodes[colnames(original_mat)[3:ncol(original_mat)]]
  
  #subsample to make it quicker for testing (remove later)
  #original_mat <- original_mat[sample(rownames(original_mat), size = 1000),]
  
  #now generate an output matrix for each strain as well as one with all cells
  for(Strain in c("HU3", "BPK081", "both")){
    
    if(Strain != "both"){
      cols <- cells_meta %>%
        filter(sample == mat_sample, strain == Strain, cell %in% colnames(original_mat)) %>%
        pull("cell")
      mat <- original_mat[,c("chromosome", "position", cols)]
    }else{
      mat <- original_mat
    }
    
    #remove kDNA
    mat <- mat %>%
      filter(chromosome != "Ld37")
    
    #convert the matrix to a proportion matrix
    prop_mat <- mat
    prop_mat[, 3:ncol(prop_mat)] <- lapply(prop_mat[, 3:ncol(prop_mat)], function(col) {
      sp <- strsplit(col, ",", fixed = TRUE)
      ref <- vapply(sp, function(z) as.integer(z[[1]]), integer(1L))
      alt <- vapply(sp, function(z) as.integer(z[[2]]), integer(1L))
      af <- alt / (ref + alt)
      return(af)
    }) 
    
    #remove unvariable loci
    unvar_loci <- rowMeans(prop_mat[,3:ncol(prop_mat)], na.rm = TRUE)
    unvar_loci <- which(is.na(unvar_loci) | unvar_loci %in% c(0,1))
    mat <- mat[-unvar_loci,]
    prop_mat <- prop_mat[-unvar_loci,]
    
    #remove loci where alternative allele is found in only 1 cells
    loci <- rowSums(prop_mat[,3:ncol(prop_mat)] > 0, na.rm = TRUE) 
    loci <- which(loci > 1)
    mat <- mat[loci,]
    prop_mat <- prop_mat[loci,]
    
    #replace heterozygous loci with reference
    mat[, 3:ncol(mat)] <- lapply(mat[, 3:ncol(prop_mat)], function(col) {
      sp <- strsplit(col, ",", fixed = TRUE)
      ref <- vapply(sp, function(z) as.integer(z[[1]]), integer(1L)) #takes left character
      alt <- vapply(sp, function(z) as.integer(z[[2]]), integer(1L)) #takes right character
      
      #convert heterozygous loci to homozygous reference
      ref <- ifelse(alt > 0 & ref > 0, alt + ref, ref) 
      alt <- ifelse(alt > 0 & ref > 0, 0, alt)
      
      #if loci with homozygous alternative has only one read, assume that is a mistake and replace by reference instead
      alt_2 <- ifelse(alt == 1 & ref == 0, 0, alt)
      ref <- ifelse(alt == 1 & ref == 0, 1, ref)
      alt <- alt_2
      rm(alt_2)
      #paste it back
      col <- paste(ref, alt, sep = ",")
      return(col)
    }) 
    
    #remove rows with too many NA values
    NA_to_remove <- data.frame(row = rownames(prop_mat),
                               position = prop_mat$position, 
                               chromosome = prop_mat$chromosome, 
                               prop_NA = rowMeans(is.na(prop_mat[,c(3:ncol(prop_mat))]))) %>%
      mutate(to_remove = prop_NA >= NA_threshold)
    
    #make a plot showing the bins removed due to excessive NA values
    plot <- NA_to_remove %>%
      ggplot(aes(y = prop_NA, x = position, color = to_remove))+
      geom_point(size = 0.1)+
      ggtitle(Strain)+
      facet_wrap(vars(chromosome), scales = "free_x", ncol = 1)+
      theme_void()+
      theme(axis.text.y = element_text(), axis.title.y = element_text())
    
    dens_plot <- NA_to_remove %>%
      #filter(chromosome %in% "Ld01") %>%
      ggplot(aes(y = prop_NA, x = after_stat(density))) +
      stat_density(orientation = "y", geom = "area", position = "identity", alpha = 0.6, fill = "lightblue") +
      #stat_density(aes(fill = to_remove), orientation = "y", geom = "area", position = "identity", alpha = 0.6) +
      geom_hline(yintercept = NA_threshold, linetype = "dashed")+
      facet_wrap(vars(chromosome), ncol = 1, scales = "free_y") +
      theme_void() + guides(fill = "none")
    
    plot <- ggalign::align_plots(plot, dens_plot, nrow = 1, guides = "r")
    
    #remove positions with excessive NA values
    rows <- NA_to_remove %>%
      filter(prop_NA <= NA_threshold) %>%
      pull("row")
    
    mat <- mat[rows,]
    
    #convert it to a matrix-like dataframe
    mat <- mat %>%
      as.data.frame() %>%
      mutate(row_name = paste0(chromosome, "_", position)) %>%
      column_to_rownames("row_name") %>%
      select(-chromosome, -position) 
    
    #export the plot
    ggsave(plot, file = paste0("inputs/nucleotide_variants/scistree2_data/",mat_sample, "_", Strain, "_NA_removal_plot.pdf"), height = 60, width = 5, units = "in", limitsize = FALSE)
    
    #export the matrix
    file_name <- paste0(mat_sample, "_", Strain, "_ad_matrix_for_scistree2.tsv")
    mat %>%
      rownames_to_column(var = "site") %>%
      write_tsv(file = paste0("inputs/nucleotide_variants/scistree2_data/", file_name), col_names = TRUE) 
  }
}


# PARALLEL over files
furrr::future_map(matrices, process_one_matrix, .progress = TRUE)