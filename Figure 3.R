#clean the environment####
rm(list = ls(envir = parent.frame()))
gc()

#load libraries####
library(tidyverse)
library(ggridges)
library(ggrepel)
library(patchwork)
library(circlize)
library(ComplexHeatmap)
library(tidyHeatmap)
source("scDNA_functions_S3.R")

#get inputs####
true_cells <- readRDS("inputs/true_cells.rds")
scDNA10X <- readRDS("inputs/scDNA10X.rds")

##set base theme parameters for ggplot####
update_geom_defaults("point", list(size = 0.5))

#bind cells meta####
cells_meta <- lapply(c(true_cells, scDNA10X), function(x)x$metadata$cells_meta)
cells_meta <- lapply(cells_meta, function(x){
  x$sample <- as.character(x$sample) 
  return(x)})
cols = Reduce(intersect, lapply(cells_meta, names))
cells_meta <- lapply(cells_meta, function(x)x[,cols])
cells_meta <- bind_rows(cells_meta)

#plot the figures####
##plot figure 2A####
figure_2A <- cells_meta %>%
  mutate(sample = ifelse(sample == "10X", paste(sample, strain), paste("Sample", sample))) %>%
  select(sample, strain, ICCV, ICF_score, mean_int_dist, max_int_dist, effective_depth_of_coverage) %>%
  pivot_longer(cols = -c(sample, strain), names_to = "metric") %>%
  filter(metric != "max_int_dist") %>%
  mutate(sample = factor(sample, levels = c("Sample 1","Sample 2", "Sample 3", "Sample 5", "Sample 6", "10X BPK081", "10X HU3"))) %>%
  ggplot(aes(y = fct_rev(sample), x = value, fill = sample))+
  geom_density_ridges(alpha = 0.6, scale = 2)+
  labs(x = NULL, y = NULL)+
  guides(fill = "none")+
  facet_wrap(vars(metric), scales = "free_x", ncol = 1)

##plot_figure 2B####
###creat a color pallete
colors <- colorRamp2(breaks = c(0:10), colors = heat_col(c(0:10)))

#create the plot list where both plots will be stored
plot_list <- list()
#make the same plot for samples 1 (index 1) and sample 6 (index 5 because sample 4 was removed from true_cells)
for(i in c(1, 5)){
  #remove outlier bins
  bins <- true_cells[[i]]$metadata$bins_meta %>%
    filter(!is_outlier & is_mappable) %>%
    rownames()
  
  #get normalized counts
  matrix_to_plot <- true_cells[[i]]$counts$corrected_counts[bins,]
  matrix_to_plot <- apply(matrix_to_plot, 2, normalize, method = "density_peak")
  
  #prepare annotations
  sample <- as.character(unique(true_cells[[i]]$metadata$cells_meta$sample))
  annotation_colors <- list(strain = c(BPK081 = "#b7776e", HU3 = "#6eb777", doublet = "#776eb7", low_cor = "lightgrey"),
                            is_outlier = c(`TRUE` = "orange", `FALSE` = "lightgrey"))
  cell_annotation <- HeatmapAnnotation(df = true_cells[[i]]$metadata$cells_meta[colnames(matrix_to_plot),c("strain", "is_outlier")], 
                                       col = annotation_colors, show_legend = i != 1, annotation_name_gp = gpar(fontsize = 10))
  
  #breake rows per chromosome
  row_breaks <- true_cells[[i]]$metadata$bins_meta %>%
    rownames_to_column("bin_id") %>%
    select(bin_id, chromosome) %>%
    deframe()
  row_breaks <- row_breaks[rownames(matrix_to_plot)]
  
  plot <- Heatmap(matrix_to_plot*2, 
                  show_row_names = FALSE,
                  show_column_names = FALSE,
                  cluster_rows = FALSE,
                  col = colors, 
                  row_split = row_breaks,
                  row_title_rot = 0,
                  row_title_gp = gpar(fontsize = 8),
                  show_heatmap_legend = i != 1,
                  name = "count",
                  top_annotation = cell_annotation,
                  use_raster = TRUE,
                  gap = unit(0.5, "mm"),
                  raster_quality = 15,
                  column_title = paste("Sample", sample),
                  clustering_method_columns = "ward.D2")
  
  plot_list <- c(plot_list, plot)
}

plot1 <- wrap_elements(full = grid.grabExpr(draw(plot_list[[1]])))
plot2 <- wrap_elements(full = grid.grabExpr(draw(plot_list[[2]])))

figure_2B <- plot1 + plot2

final <- figure_2A + figure_2B + 
  plot_layout(widths = c(0.15, 0.85))+
  plot_annotation(tag_levels = 'A')&
  theme(plot.tag = element_text(size = 18, face = "bold"),
        plot.tag.position   = c(0, 1),    # top-left corner
        plot.tag.background = element_rect(fill = NA, colour = NA),
        axis.title = element_text(size = 10),
        text = element_text(size = 10))

#save the final figure
plot_scale <- 1.8
ggsave("figure_2.pdf", plot = final, width = 8.27 * plot_scale, height = 5 * plot_scale)

#open it
system2('open', args = "figure_2.pdf", wait = FALSE)