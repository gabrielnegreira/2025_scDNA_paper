#clean the environment####
rm(list = ls(envir = parent.frame()))
gc()


#this script will read the raw allele depth matrices, clean them and export them so they can be used with the `ScisTree2.py` script.
#libraries
library(tidyverse)
library(readr)
library(UpSetR)
library(furrr)

# pick a sensible number of workers
plan(multisession, workers = parallel::detectCores() - 1)

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

#create a function to process the matrices so it can run in parallel
process_one_matrix <- function(mfile) {
  #get matrix
  original_mat <- readr::read_tsv(file.path("inputs/nucleotide_variants", mfile), show_col_types = FALSE)
  
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
    
    #convert the matrix to a proportion matrix
    prop_mat <- mat
    prop_mat[, 3:ncol(prop_mat)] <- lapply(prop_mat[, 3:ncol(prop_mat)], function(col) {
      sp <- strsplit(col, ",", fixed = TRUE)
      ref <- vapply(sp, function(z) as.integer(z[[1]]), integer(1L))
      alt <- vapply(sp, function(z) as.integer(z[[2]]), integer(1L))
      alt / (ref + alt)
    }) 
    
    #remove unvariable loci
    unvar_loci <- rowMeans(prop_mat[,3:ncol(prop_mat)], na.rm = TRUE)
    unvar_loci <- which(is.na(unvar_loci) | unvar_loci %in% c(0,1))
    mat <- mat[-unvar_loci,]
    
    #convert it to a matrix-like dataframe
    r_names <- paste(mat$chromosome, mat$position, sep = "_")
    mat <- mat[,c(3:ncol(mat))]
    mat <- as.data.frame(mat)
    rownames(mat) <- r_names
    rm(r_names)
    
    #export the matrix
    file_name <- paste0(mat_sample, "_", Strain, "_ad_matrix_for_scistree2.tsv")
    mat %>%
      rownames_to_column(var = "site") %>%
      write_tsv(file = paste0("inputs/nucleotide_variants/", file_name), col_names = TRUE) 
  }
}


# PARALLEL over files
furrr::future_map(matrices, process_one_matrix, .progress = TRUE)