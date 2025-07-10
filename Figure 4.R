#clean the environment####
rm(list = ls(envir = parent.frame()))
gc()

#load libraries####
source("color_palettes.R")
library(tidyverse)
library(ggridges)
library(ggalign)

#get inputs
true_cells <- readRDS("inputs/true_cells_1kb_bins.rds")
cnv_files <- list.files("inputs", pattern = "100nt_CNV")
cnv_files <- lapply(cnv_files, function(x)read.csv(paste0("inputs/", x), row.names = 1))
cnv_files <- lapply(cnv_files, function(x){
  colnames(x) <- gsub("_", "", colnames(x))
  return(x)
})

#plot
heatmaps <- list()
densities <- list()
for(i in seq_along(cnv_files)){
   cnv_mat <- cnv_files[[i]]
   cnv_name <- names(cnv_files)[i]
   sample <- 6
   index <- ifelse(sample >= 4, sample - 1, sample)
   strains <- true_cells[[index]]$metadata$cells_meta %>%
     rownames_to_column("cell") %>%
     select(cell, strain) %>%
     deframe()
   names(strains) <- gsub("_.*", "", names(strains))
   strains <- strains[colnames(cnv)]
   
   
   cnv[,] <- apply(cnv, 2, normalize, method = "lipari")
   ggheatmap(cnv)+
    scale_fill_gradientn(colors = rev(brewer.pal(5, name = "RdYlBu")))+
    #scale_fill_gradientn(colors = scico(30, palette = "imola"))+
    scale_x_discrete(breaks = NULL)+
    anno_top(size = 0.2)+
    align_dendro()+
    anno_top() + #create a second annotation space, below the dendogram and on top of the heatmap
    ggalign(data = strains[colnames(cnv)], size = 0.1)+ #initiate a ggplot object withy the strains as data.
    geom_tile(aes(y = 1, fill = factor(value))) + #add annotation tiles.
    scale_fill_manual(name = "Strain", values = strain_colors)+ #set the annotation colors
    theme_void()
}