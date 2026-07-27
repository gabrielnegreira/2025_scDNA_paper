#clean the environment####
source("clean_environment.R")

#load libraries####
library(tidyverse)
library(ggridges)
library(ggalign)
library(ggpubr)
library(rstatix)
source("git_modules/scDNA_tools/scDNA_functions.R")
source("color_palettes.R")

#get inputs####
all_SPCs <- readRDS("inputs/karyotyping_objects/all_SPCs.rds")
true_cells <- readRDS("inputs/karyotyping_objects/true_cells.rds")
scDNA10X <- readRDS("inputs/karyotyping_objects/scDNA10X.rds")
cells_meta <- read_delim("inputs/cell_qc/cells_meta.tsv") %>%
  column_to_rownames("rowname") %>%
  filter(cell_or_background == "cell" & strain != "doublet" & sample != "SPC-STD4")
all_SPCs[["10X"]] <- scDNA10X
true_cells[["10X"]] <- scDNA10X
rm(scDNA10X)
##set base theme parameters for ggplot####
update_geom_defaults("point", list(size = 0.1))
update_geom_defaults("boxplot", list(outlier.size = 0.5))

#fix sample names
names(all_SPCs) <- sample_names[names(all_SPCs)]
names(true_cells) <- sample_names[names(true_cells)]

#bind bins_meta####
bins_meta <- lapply(true_cells, function(x) x$metadata$bins_meta)
bins_meta <- bind_rows(bins_meta) 

bins_meta <- bins_meta %>%
  mutate(sample = ifelse(experiment == "10X", paste(sample, strain), sample)) %>%
  mutate(sample = factor(sample_names[sample], levels = sample_names))

cells_meta <- cells_meta %>%
  mutate(sample = ifelse(experiment == "10X", paste(sample, strain), sample)) %>%
  mutate(sample = factor(sample_names[sample], levels = sample_names))

#plot the figures####
#create an empty list to store the figures
figures <- list()

##plot figure 2A####
figures[[length(figures) + 1]] <- bins_meta %>%
  #filter(experiment == "Atrandi") %>%
  select(sample, gc_content, mean_raw_counts, mean_corrected_counts) %>%
  filter(sample != "SPC-STD4") %>%
  pivot_longer(cols = c(mean_raw_counts, mean_corrected_counts), values_to = "count", names_to = "count_type") %>%
  mutate(count_type = c(mean_raw_counts = "Raw Counts", mean_corrected_counts = "Corrected Counts")[count_type]) %>%
  group_by(sample, count_type) %>%
  filter(!count %in% boxplot.stats(count)$out) %>%
  ggplot(aes(x = gc_content, y = count))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE, color = "blue", linewidth = 1) +  # Linear model line
  scale_x_continuous(limits = c(0.5, 0.7))+
  labs(x = "GC content", y = "Read Counts")+
  facet_grid(rows = vars(sample), cols = vars(fct_rev(count_type)), scales = "free")+
  theme_bw()+
  theme(panel.grid.minor = element_blank())

##plot figure 2B####
#compare raw vs corrected counts
to_plot <- list()
to_plot[["raw"]] <- lapply(true_cells, function(x)x$counts$raw_counts)
breaks <- unique(bins_meta$chromosome)
breaks <- breaks[seq(1, 35, by = 2)]

figures[[length(figures) + 1]] <- bins_meta %>%
  #filter(experiment == "Atrandi") %>%
  group_by(sample) %>%
  mutate(bin_position = row_number()) %>%
  filter(!is_outlier & !is_empty) %>%
  filter(sample != "SPC-STD4") %>%
  pivot_longer(cols = c("mean_raw_counts", "mean_corrected_counts"), names_to = "count_type", values_to = "count") %>%
  group_by(sample, chromosome) %>%
  filter(!count %in% boxplot.stats(count)) %>%
  mutate(count_type = c(mean_raw_counts = "Raw Counts", mean_corrected_counts = "Corrected Counts")[count_type]) %>%
  ggplot(aes(x = chromosome, y = count))+
  geom_boxplot(outliers = FALSE)+
  scale_x_discrete(breaks = breaks, labels = seq(1, 35, by = 2))+
  guides(color = "none")+
  labs(x = "Chromosome", y = "Mean count")+
  facet_grid(cols = vars(fct_rev(count_type)), rows = vars(sample), scale = "free")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())

##plot figure 2C#####
##compare lorenz curves with and without GC correction
###since the 10X data has more rows than the Atrandi data, will group rows into 100 bins and take their average.
n_bins <-1000
lorenz_df <- lapply(true_cells, function(x){
  #start by creating a list with the raw and corrected counts for each dataset
  raw_counts <- x$counts$raw_counts
  corrected_counts <- x$counts$corrected_counts
  count_list <- list(raw = raw_counts, corrected = corrected_counts)
  
  #now for each count matrix...
  count_list <- lapply(count_list, function(mat){
    #first remove rows which have no counts (specially in 10X dataset, where artificial Ns were introduced in reference genome)
    not_empty <- rowSums(mat)
    not_empty <- not_empty[not_empty > 0]
    mat <- mat[names(not_empty),]
    
    #then group rows into 1000 bins...
    n_rows <- nrow(mat)
    bins <- cut(seq_len(n_rows), breaks = n_bins, labels = FALSE)
    names(bins) <- 1:n_rows
    new_mat <- matrix(nrow = n_bins, ncol = ncol(mat))
    colnames(new_mat) <- colnames(mat)
    #then calcula the means in each row group
    for(row_group in unique(bins)){
      rows <- as.integer(names(bins[bins == row_group]))
      new_mat[row_group,] <- colMeans(mat[rows,,drop = FALSE])
    }
    #now for each cell (col), calculate the sorted cumulative sum
    new_mat <- apply(mat, 2, function(col){
      col <- sort(col)
      col <- cumsum(col)
      col <- col/max(col)
      return(col)
    })
    #normalize the matrix
    #new_mat[,] <- apply(new_mat, 2, function(x)x/sum(x))
    
    #convert the matrix to long data frame, compute the mean count and the standard error
    new_mat <- as.data.frame(new_mat) %>%
      rownames_to_column("frac_genome") %>%
      mutate(frac_genome = as.integer(frac_genome)) %>%      
      pivot_longer(cols = -frac_genome, names_to = "cell", values_to = "frac_coverage") %>%
      group_by(frac_genome) %>%
      summarise(mean = mean(frac_coverage), sd = sd(frac_coverage), ncells = n()) %>%
      mutate(se = sd / sqrt(ncells)) %>%
      mutate(frac_genome = frac_genome/max(frac_genome))
    return(new_mat)
  })
  #add to the data frame which data type it is
  count_list$raw$data <- "raw"
  count_list$corrected$data <- "corrected"
  x <- bind_rows(count_list)
  return(x)
})

#add to each data frame, from which sample it originates.
for(Sample in names(lorenz_df)){
  lorenz_df[[Sample]]$sample <- Sample
}

#save everything as a single data frame (for plotting)
lorenz_df <- bind_rows(lorenz_df) %>%
  mutate(experiment = ifelse(grepl("10X", sample), "10X", "Atrandi"))

#calculate inflection point (for placing labels in the plot)
labels_df <- lorenz_df %>%
  group_by(data, sample) %>%
  summarise(x = find_knee(x = frac_genome, y = mean, val_to_return = "x"),
            y = find_knee(x = frac_genome, y = mean, val_to_return = "y")) %>%
  filter(data == "corrected") %>%
  filter(sample %in% c("10X", "SPC-STD2", "SPC-PTA2")) %>%
  mutate(label = ifelse(sample == "SPC-STD2", "SPC-STD1-3",
                        ifelse(sample == "SPC-PTA2", "SPC-PTA1-2",
                               sample)))
  
figures[[length(figures) + 1]] <- lorenz_df %>%
  mutate(data = fct_rev(factor(data))) %>%
  filter(data == "corrected") %>%
  #filter(experiment == "Atrandi") %>%
  #mutate(mean = frac_coverage) %>%
  ggplot(aes(x = frac_genome, y = mean, color = sample))+
  geom_line(size = 0.2)+
  geom_line(data = data.frame(frac_genome = c(0,1), mean = c(0,1)), color = "lightgrey", linetype = "dashed")+
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se, group = sample), alpha = 0.2, color = NA)+
  geom_text_repel(data = labels_df, aes(x = x, y = y, label = label), show.legend = FALSE, force = 30, nudge_x = 0.2, size = 2, min.segment.length = 0)+
  #facet_wrap(vars(sample), ncol = 1)+
  guides(color = "none")+
  labs(x = "Fraction of the genome", y = "cumulative fraction of counts")+
  scale_color_manual(values = sample_colors)+
  theme_bw()+
  theme(panel.grid = element_blank())

##plot figure 2D####
plot_list <- list()
min_ncells <- cells_meta %>%
  mutate(sample = ifelse(experiment == "10X", "10X", sample)) %>%
  group_by(sample) %>%
  summarise(n_cells = n()) %>%
  pull("n_cells") %>%
  min()

for(metric in c("gini", "mapd", "ICF_score")){
  stat_data <- cells_meta %>%
    mutate(sample = ifelse(experiment == "10X", "10X", as.character(sample))) #%>%
    #slice_sample(n = min_ncells, by = sample)
  
  stat_test <- stat_data %>%
    mutate(sample = as.character(sample)) %>%
    wilcox_test(formula = as.formula(paste0(metric, " ~ sample")), ref.group = "10X", p.adjust.method = "hochberg") %>%
    add_significance("p.adj")

  y_positions <- stat_data %>%
    group_by(sample) %>%
    summarise(y.position = max(.data[[metric]], na.rm = TRUE)) %>%
    rename(group2 = sample)
  
  stat_test <- stat_test %>%
    left_join(y_positions, by = "group2") %>%
    mutate(y.position = y.position * 1.1)
  
  plot_list[[metric]] <- cells_meta %>%
    mutate(plot_label = metric) %>%
    ggplot(aes(x = factor(sample), y = .data[[metric]]))+
    geom_violin(aes(fill = factor(sample)), size = 0.2)+
    geom_boxplot(width = 0.075, outlier.size = 0.1, size = 0.2)+
    geom_text(data = stat_test, aes(x = group2, y = y.position, label = p.adj.signif))+
    #stat_pwc(aes(group = sample), ref.group = "10X", method = "wilcox_test", label = "p.adj.signif", p.adjust.method = "hochberg", hide.ns = TRUE)+
    #stat_pvalue_manual(stat_test, label = "p.adj.signif")+
    labs(x = NULL, y = c(gini = "Gini", mapd = "MAPD", "ICF_score" = "S-phase Score")[metric])+
    guides(fill = "none")+
    scale_y_continuous(expand =c(0, 0.15))+
    scale_fill_manual(values = sample_colors)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
}

#combine them
figures[[length(figures) + 1]] <- ggalign::align_plots(!!!plot_list, nrow = 1)

##plot_figure 2E####
###creat a color pallete
breaks <- lapply(true_cells, function(x)x$counts$normalized_counts) #get the normalized counts of each sample
breaks <- breaks[c(1,5)] #subset only samples 1 and 5
breaks <- unique(as.integer(unlist(breaks)))
breaks <- c(0:max(breaks))
color_palette <- heat_col(breaks)
names(color_palette) <- breaks
#create the plot list where both plots will be stored
plot_list <- list()
#make the same plot for samples 2 (index 2) and SPC-PTA2 (index 5 because SPC-STD4 was removed from true_cells)
for(i in c("10X", "SPC-STD2", "SPC-PTA2")){
  #remove outlier bins
  bins <- true_cells[[i]]$metadata$bins_meta %>%
    filter(!is_outlier & is_mappable) %>%
    rownames_to_column("bin_id") %>%
    select(bin_id, chromosome) %>%
    deframe() %>%
    factor() %>%
    fct_rev()
  
  #remove outlier cells
  cells <- true_cells[[i]]$metadata$cells_meta %>%
    filter(!is_outlier) %>%
    rownames()
  
  #get normalized counts
  matrix_to_plot <- true_cells[[i]]$counts$normalized_counts[names(bins), cells] #get corrected counts
  matrix_to_plot <- apply(matrix_to_plot, 2, normalize, method = "density_peak")
  matrix_to_plot <- matrix_to_plot*2
  #get cells' strain
  strains <- true_cells[[i]]$metadata$cells_meta %>%
    rownames_to_column("cell") %>%
    select(cell, strain) %>%
    deframe()
  
  #plot density
  dens_plot <- ggplot(data.frame(value = as.vector(matrix_to_plot)), aes(x = value))+
    geom_density()+
    scale_x_continuous(name = "Normalized Read Count\n(multiplied by 2)", breaks = c(1:100), limits = c(1,12))
  
  #plot core heatmap
  hm_plot <- ggheatmap(matrix_to_plot)+#, filling = NULL) + #start the heatmap
    #geom_tile(aes(fill = value))+
    #raster_magick(geom_tile(aes(fill = value)), res = 1200) +
    scale_fill_gradientn(name = "Normalized Read Count\n(multiplied by 2)", colors = color_palette, breaks = breaks, limits = c(0, max(breaks))) + #set the fill scale
    scale_x_discrete(labels = NULL)+ #remmoves cell labels
    scale_y_discrete(labels = NULL, breaks = NULL)+ #removes bin labels
    facet_grid(switch = "y")+ #this will place the row group labels on the left side instead of the right side. 
    theme(strip.text.y.left = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 6), #this is needed to keep panel labels (chromosomes) as default `theme_ggalign()` removes them. 
          strip.placement = "outside", #allow labels to be draw outside the box
          strip.clip = "off")+ #prevent labels from being clipped. 
    anno_left(size = 0) + #start an empty anotation box to the left (used to group rows per chromosome)
    align_group(group = bins[rownames(matrix_to_plot)])+ #group bins by chromosome
    anno_top(size = 0.3) + #create an anotation space in the top
    free_border(dens_plot, borders = "l") + #initialize a ggplot object for the density plot
    patch_titles(top = i)+ #add a title to the top annotation
    anno_top() + #create another annotation space at the top
    align_dendro(method = "ward.D2", size = 0.2)+ #add the dendogram (which also reorders the columns)
    theme_void()+ #remove irrelevant elements of the dendogram plot
    anno_top() + #create a second annotation space, below the dendogram and on top of the heatmap
    ggalign(data = strains[colnames(matrix_to_plot)], size = 0.2)+ #initiate a ggplot object withy the strains as data.
    geom_tile(aes(y = 1, fill = factor(value))) + #add annotation tiles.
    scale_fill_manual(name = "Strain", values = strain_colors)+ #set the annotation colors
    theme_void() & #remove irrelevant elements of the annotation box
    theme(panel.spacing = unit(1, "pt"), 
          legend.key.size = unit(15, "pt"),
          legend.position = "bottom")  
  plot_list <- c(plot_list, hm_plot) 
}

#align both figures
figures[[length(figures) + 1]] <- ggalign::align_plots(!!!plot_list, guides = "rltb")

top <- ggalign::align_plots(!!!figures[c(1, 2)], widths = c(0.2, 0.8))
middle <- ggalign::align_plots(free_border(figures[[3]], borders = "b"), figures[[4]], widths = c(0.2, 0.8))
final <- ggalign::align_plots(middle, figures[[5]], ncol = 1, heights = c(0.25, 0.75))
final <- final + layout_tags("A") + layout_theme(plot.tag = element_text(size = 16))

#save the final figure panel####
plot_scale <- 1.8
ggsave("figure_2.pdf", plot = final, width = 8.27 * plot_scale, height = 7 * plot_scale)

#open it
if (Sys.info()["sysname"] == "Darwin") {
  system2("open", args = "figure_2.pdf", wait = FALSE)
} else if (Sys.info()["sysname"] == "Linux") {
  system2("xdg-open", args = "figure_2.pdf", wait = FALSE)
}
