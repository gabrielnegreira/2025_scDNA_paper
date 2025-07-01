#clean the environment####
rm(list = ls(envir = parent.frame()))
gc()

#load libraries####
library(tidyverse)
library(patchwork)

#get inputs####
true_cells <- read_rds("inputs/true_cells.rds")
plot_list <-list()
for(i in c(1:5)){
  cells_meta <- true_cells[[i]]$metadata$cells_meta
  rownames(cells_meta) <- gsub("_.*", "", rownames(cells_meta))
  
  sample <- ifelse(i >= 4, i + 1, i) #since sample 4 is removed, indexes changed. 
  af_matrix <- read_tsv(paste0("inputs/sample_", sample, "_af_matrix_formated.tsv")) %>%
    mutate(row_name = paste0(chromosome, "_", position)) %>%
    select(-chromosome, -position) %>%
    column_to_rownames("row_name") %>%
    as.matrix()
  
  #imput NA values
  pca_mat <- apply(af_matrix, 1, function(x)ifelse(is.na(x), mean(x, na.rm = TRUE), x))
  
  #try PCA first
  pca <- prcomp(pca_mat)
  plot_list[[i]] <- cbind(pca$x, cells_meta[rownames(pca$x), ]) %>%
   ggplot(aes(x = PC1, y = PC2, color = strain))+
   geom_point()+
   ggtitle(paste("Sample", sample))
  rm(pca_mat)

}

wrap_plots(plot_list)

#now do it per strain
plot_list <-list()
for(i in c(1:5)){
  cells_meta <- true_cells[[i]]$metadata$cells_meta
  rownames(cells_meta) <- gsub("_.*", "", rownames(cells_meta))
  sample <- ifelse(i >= 4, i + 1, i) #since sample 4 is removed, indexes changed. 
  af_matrix <- read_tsv(paste0("inputs/sample_", sample, "_af_matrix_formated.tsv")) %>%
    mutate(row_name = paste0(chromosome, "_", position)) %>%
    select(-chromosome, -position) %>%
    column_to_rownames("row_name") %>%
    as.matrix()
  
  plot_each_strain <- list()
  for(Strain in c("BPK081", "HU3")){
    #imput NA values
    cells <- cells_meta %>%
      filter(strain == Strain) %>%
      rownames() %>%
      gsub("_.*", "", .)
    
    #get matrix for strain cells only
    pca_mat <- af_matrix[,cells]
    
    #imput NAs
    pca_mat <- t(apply(pca_mat, 1, function(x)ifelse(is.na(x), mean(x, na.rm = TRUE), x)))
    
    #remove NAs-only rows
    rows_to_keep <- rowSums(is.na(pca_mat))
    rows_to_keep <- names(rows_to_keep[rows_to_keep < ncol(pca_mat)])
    pca_mat <- pca_mat[rows_to_keep,]
    
    #try PCA first
    pca <- prcomp(t(pca_mat))
    plot_each_strain[[Strain]] <- cbind(pca$x, cells_meta[rownames(pca$x), ]) %>%
      ggplot(aes(x = PC1, y = PC2, color = strain))+
      geom_point()+
      ggtitle(paste0("Sample: ", sample, "; Strain: ", Strain))
    rm(pca_mat)
  }
  plot_list[[i]] <- wrap_plots(plot_each_strain, nrow = 1)
}

wrap_plots(plot_list, ncol = 1)


