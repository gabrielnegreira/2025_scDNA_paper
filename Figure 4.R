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
for(cnv in cnv_files){
 ggheatmap(cnv)+
    scale_x_discrete(breaks = NULL)+
    anno_top(size = 0.2)+
    align_dendro()+
    anno_top() + #create a second annotation space, below the dendogram and on top of the heatmap
    ggalign(data = strains[colnames(matrix_to_plot)], size = 0.1)+ #initiate a ggplot object withy the strains as data.
    geom_tile(aes(y = 1, fill = factor(value))) + #add annotation tiles.
    scale_fill_manual(name = "Strain", values = strain_colors)+ #set the annotation colors
    theme_void()
}