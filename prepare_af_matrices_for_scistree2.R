#clean the environment####
rm(list = ls(envir = parent.frame()))
gc()


#this script will read the raw allele depth matrices, clean them and export them so they can be used with the `ScisTree2.py` script.
#libraries
library(tidyverse)
library(readr)
library(UpSetR)

#get the data
matrices <- list.files("inputs/nucleotide_variants/", pattern = "ad_matrix.tsv")

#get cells_metadata (needed for spliting matrices by strain)
cells_meta <- read_delim("inputs/cell_qc/cells_meta.tsv") %>%
  column_to_rownames("rowname") %>%
  rownames_to_column("cell") %>%
  select(sample, strain, cell) %>%
  mutate(sample = tolower(sample)) %>%
  mutate(sample = gsub(" ", "_", sample)) %>%
  mutate(simple_barcode = gsub("_.*", "", cell))


for(i in seq_along(matrices)){
  #get matrix
  original_mat <- read_tsv(paste0("inputs/nucleotide_variants/", matrices[i]))

  #remove last column containing only NA values (issue with tsv file)
  cols_to_remove <- apply(original_mat, 2, function(x) all(is.na(x)))
  cols_to_keep <- which(!cols_to_remove)
  original_mat <- original_mat[,cols_to_keep]
  
  #specify to which sample this matrix comes from 
  mat_sample <- gsub("_ad_matrix.*", "", matrices[i])
  
  #create a key name pair for the correct barcodes
  correct_barcodes <- cells_meta %>%
    filter(sample == mat_sample) %>%
    select(simple_barcode, cell) %>%
    deframe()
  
  #check if the colnames match the names in the cells_meta
  # upset(fromList(list(cells_meta = names(correct_barcodes),
  #                     colnames = colnames(original_mat[,3:ncol(original_mat)]))))
  
  #keep only cells found in the cells_meta
  original_mat <- original_mat[,which(colnames(original_mat) %in% c("chromosome", "position", names(correct_barcodes)))]
  
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
    
    #convert the matrix to a count matrix
    count_mat <- mat
    count_mat[,3:ncol(mat)] <- apply(count_mat[,3:ncol(count_mat)], 2, function(x){
      x <- sapply(x, function(y){
        y <- strsplit(y, split = ",")
        y <- unlist(y)
        y <- as.integer(y)
        y <- sum(y, na.rm = TRUE)
        return(y)
      })
    })
    
    #now count how many counts each bin had and save it in a dataframe
    af_counts <- data.frame(chromosome = count_mat$chromosome, 
                            position = count_mat$position, 
                            total_counts = (rowSums(count_mat[,3:ncol(count_mat)])), 
                            row = rownames(count_mat))
    rm(count_mat)
    #get the top 1000 loci for each chromosome
    rows <- af_counts %>%
      group_by(chromosome) %>%
      slice_max(n = 1000, order_by = total_counts) %>%
      pull("row")
    
    #subset the matrix
    mat <- mat[rows,]
    rm(rows)
    
    #convert it to a matrix
    r_names <- paste(mat$chromosome, mat$position, sep = "_")
    mat <- mat[,c(3:ncol(mat))]
    mat <- as.data.frame(mat)
    rownames(mat) <- r_names
    rm(r_names)
    
    #export the matrix
    file_name <- paste0(mat_sample, "_", Strain, "_ad_matrix_for_scistree2.tsv")
    write_tsv(mat, file = paste0("inputs/nucleotide_variants/", file_name)) 
  }
}