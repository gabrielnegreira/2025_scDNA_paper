#clean the environment####
rm(list = ls(envir = parent.frame()))
gc()


#this script will read the raw allele depth matrices, clean them and expor them so they can be used with the `ScisTree2.py` script.
#libraries
library(tidyverse)

#get the data
matrices <- list.files("inputs/nucleotide_variants/", pattern = "ad_matrix_formated.tsv")

for(i in seq_along(matrices)){
  mat <- read_tsv(paste0("inputs/nucleotide_variants/", matrices[i]))
  #remove last column containing only NA values (issue with tsv file)
  cols_to_remove <- apply(mat, 2, function(x) all(is.na(x)))
  cols_to_keep <- which(!cols_to_remove)
  mat <- mat[,cols_to_keep]

  
  #subsample to make it quicker for testing (remove later)
  mat <- mat[sample(rownames(mat), size = 1000),]
  
  # #convert the matrix to a count matrix
  # count_mat <- mat[,3:ncol(mat)]
  # count_mat[,] <- apply(count_mat, 2, function(x){
  #   x <- sapply(x, function(y){
  #     y <- strsplit(y, split = ",")
  #     y <- unlist(y)
  #     y <- as.integer(y)
  #     y <- sum(y, na.rm = TRUE)
  #     return(y)
  #   })
  # })
  # 
  # #now count how many counts each bin had
  # rowSums(count_mat)
  
  #convert it to long tible
  alleles <- mat %>%
    rownames_to_column("rowname") %>%
    pivot_longer(cols = -c(chromosome, position, rowname), names_to = "cell", values_to = "allele_depth") %>%
    separate(allele_depth, into = c("ref_count", "alt_count"), sep = ",", remove = FALSE, convert = TRUE) %>%
    mutate(cell = factor(cell)) %>% #prevents cells from being completely removed when converting the tibble back to the matrix
    mutate(total_count = ref_count + alt_count,
           allele_frequency = alt_count/total_count,
           variant_type = ifelse(allele_frequency == 0 | allele_frequency == 1, "homozygous", "heterozygous")) %>%
    filter(variant_type == "homozygous") %>% #remove heterozygous loci
    #filter(total_count >= 2) %>% #keep only loci with at least 2 counts
    group_by(chromosome, position) %>%
    mutate(loci_total_count = sum(total_count, na.rm = TRUE), 
              total_NAs = sum(is.na(total_count)), 
              ncells = n()) %>%
    group_by(chromosome)  %>%
    filter(loci_total_count <= sort(loci_total_count, decreasing = TRUE)[1000]) %>%
    select(rowname, chromosome, position, cell, allele_depth) %>%
    pivot_wider(names_from = "cell", values_from = "allele_depth")
    
}