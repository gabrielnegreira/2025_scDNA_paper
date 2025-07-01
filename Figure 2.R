#clean the environment####
rm(list = ls(envir = parent.frame()))
gc()

#load libraries####
library(tidyverse)
library(ggridges)
library(patchwork)
library(ggalign)
source("scDNA_functions_S3.R")
source("color_palettes.R")

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
  scale_fill_manual(values = sample_colors)+
  facet_wrap(vars(metric), scales = "free_x", ncol = 1)

##plot_figure 2B####
###creat a color pallete
breaks <- lapply(true_cells, function(x)x$counts$normalized_counts) #get the normalized counts of each sample
breaks <- breaks[c(1,5)] #subset only samples 1 and 5
breaks <- unique(as.integer(unlist(breaks)))
breaks <- c(0:max(breaks))
color_palette <- heat_col(breaks)
names(color_palette) <- breaks
#create the plot list where both plots will be stored
plot_list <- list()
#make the same plot for samples 1 (index 1) and sample 6 (index 5 because sample 4 was removed from true_cells)
for(i in c(1, 5)){
  #remove outlier bins
  bins <- true_cells[[i]]$metadata$bins_meta %>%
    filter(!is_outlier & is_mappable) %>%
    rownames_to_column("bin_id") %>%
    select(bin_id, chromosome) %>%
    deframe() %>%
    factor() %>%
    fct_rev()
  
  #get normalized counts
  matrix_to_plot <- true_cells[[i]]$counts$normalized_counts[names(bins),] #get corrected counts
  matrix_to_plot <- apply(matrix_to_plot, 2, normalize, method = "density_peak")
  matrix_to_plot <- matrix_to_plot*2
  #get cells' strain
  strains <- true_cells[[i]]$metadata$cells_meta %>%
    rownames_to_column("cell") %>%
    select(cell, strain) %>%
    deframe()
  
  #get the sample
  sample <- unique(true_cells[[i]]$metadata$cells_meta$sample)
  
  #plot density
  dens_plot <- ggplot(data.frame(value = as.vector(matrix_to_plot)), aes(x = value))+
    geom_density()+
    scale_x_continuous(name = "Normalized Read Count * 2", breaks = c(1:100))
  
  #plot core heatmap
  hm_plot <- ggheatmap(matrix_to_plot) + #start the heatmap
    scale_fill_gradientn(name = "Normalized Read Count * 2", colors = color_palette, breaks = breaks, limits = c(0, max(breaks))) + #set the fill scale
    scale_x_discrete(labels = NULL)+ #remmoves cell labels
    scale_y_discrete(labels = NULL, breaks = NULL)+ #removes bin labels
    facet_grid(switch = "y")+ #this will place the row group labels on the left side instead of the right side. 
    theme(strip.text.y.left = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 6), #this is needed to keep panel labels (chromosomes) as default `theme_ggalign()` removes them. 
          strip.placement = "outside", #allow labels to be draw outside the box
          strip.clip = "off")+ #prevent labels from being clipped. 
    anno_left(size = 0) + #start an empty anotation box to the left (used to group rows per chromosome)
    align_group(group = bins[rownames(matrix_to_plot)])+ #group bins by chromosome
    anno_top(size = 0.2) + #create an anotation space in the top
    free_border(dens_plot, borders = "l") + #initialize a ggplot object for the density plot
    patch_titles(top = paste("Sample", sample))+ #add a title to the top annotation
    anno_top() + #create another annotation space at the top
    align_dendro(method = "ward.D2", size = 0.2)+ #add the dendogram (which also reorders the columns)
    theme_void()+ #remove irrelevant elements of the dendogram plot
    anno_top() + #create a second annotation space, below the dendogram and on top of the heatmap
    ggalign(data = strains[colnames(matrix_to_plot)], size = 0.1)+ #initiate a ggplot object withy the strains as data.
    geom_tile(aes(y = 1, fill = factor(value))) + #add annotation tiles.
    scale_fill_manual(name = "Strain", values = strain_colors)+ #set the annotation colors
    theme_void() & #remove irrelevant elements of the annotation box
    theme(panel.spacing = unit(0.5, "pt"))  #make the space between row groups smaller.
  
  
  plot_list <- c(plot_list, hm_plot) 
}

#align both figures
plot_aligned <- ggalign::align_plots(!!!plot_list, guides = "r")

#covert it to patchwork compatible
figure_2B <- wrap_elements(full = ggalignGrob(plot_aligned))

final <- figure_2A + figure_2B + 
  plot_layout(widths = c(0.15, 0.85))+
  plot_annotation(tag_levels = 'A')&
  theme(plot.tag = element_text(size = 18, face = "bold"),
        plot.tag.position   = c(0, 1),    # top-left corner
        plot.tag.background = element_rect(fill = NA, colour = NA),
        axis.title = element_text(size = 10),
        text = element_text(size = 10))

#save the final figure panel####
plot_scale <- 1.8
ggsave("figure_2.pdf", plot = final, width = 8.27 * plot_scale, height = 5 * plot_scale)

#open it
system2('open', args = "figure_2.pdf", wait = FALSE)