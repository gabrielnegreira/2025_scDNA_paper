#clean the environment####
rm(list = ls(envir = parent.frame()))
gc()

#load libraries####
library(tidyverse)
library(ggridges)
library(ggrepel)
library(patchwork)
library(scales)
source("color_palettes.R")

#get inputs####
scDNAobj_list <- readRDS("inputs/karyotyping_objects/all_cells.rds")
cells_meta <- read_delim("inputs/cell_qc/cells_meta.tsv") %>%
  column_to_rownames("rowname")
true_cells <- readRDS("inputs/karyotyping_objects/true_cells.rds")

##set base theme parameters for ggplot####
update_geom_defaults("point", list(size = 0.5))

#bind cells meta####

#bind bins_meta####
bins_meta <- list()
for(i in c(1:length(scDNAobj_list))){
  bins_meta[[i]] <- scDNAobj_list[[i]]$metadata$bins_meta
  bins_meta[[i]]$sample <- unique(scDNAobj_list[[i]]$metadata$cells_meta$sample)
}
bins_meta <- bind_rows(bins_meta)

#create tables#####
##table summarising sequencing metrics
cells_meta %>%
  group_by(sample) %>%
  summarise(n_SPCs = n(), 
            total_count = sum(n_reads),
            median_count_per_SPC = median(n_reads))

#plot the figures####
##plot placeholder for figure 1A####

# Create an empty plot with text and border
figure_1A <- ggplot()+
  annotate("text", x = 0.5, y = 0.5, label = "Experimental Setup\nSchematic", size = 6, hjust = 0.5)+
  theme_void()+
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE)+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

##plot figure 1B####
figure_1B <- cells_meta %>%
  ggplot(aes(x = n_reads, y = fct_rev(factor(sample)), fill = factor(sample)))+
  geom_density_ridges(alpha = 0.5)+
  scale_x_continuous(transform = "log10", breaks = c(2000, 17000, 170000, 900000))+
  scale_y_discrete(name = NULL)+
  scale_fill_manual(values = sample_colors)+
  labs(x = "Total reads per cell (log10 scale)", y = "Sample")+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor = element_blank())

##plot figure 1C####
figure_1C <- cells_meta %>%
  group_by(sample) %>%
  arrange(desc(n_reads)) %>%
  mutate(rank = row_number()) %>%
  mutate(sample = paste0("Sample ", sample)) %>%
  group_by(sample, cell_or_background) %>%
  mutate(last_true_cell = rank == rank[which(rank == max(rank))]) %>%
  ungroup() %>%
  mutate(last_true_cell = last_true_cell & cell_or_background == "cell") %>%
  mutate(label = ifelse(last_true_cell, paste0("Rank: ", rank, "\nread_count: ", n_reads), NA)) %>%
  ggplot(aes(x = rank, y = n_reads, color = cell_or_background, label = label))+
  geom_line()+
  geom_text_repel(show.legend = FALSE, nudge_x = -0.5, nudge_y = -0.5, size = 2)+
  scale_x_continuous(transform = "log10")+
  scale_y_continuous(transform = "log10")+
  scale_color_manual(values = c(cell = "darkblue", backgroun = "lightgrey"))+
  guides(color = "none")+
  labs(x = "Barcode Rank (log10)", y = "Read Count (log10)")+
  facet_wrap(vars(sample), nrow = 1, scales = "free")+
  theme_minimal()

##plot figure 1D####
#check how many cells were removed from each sample
figure_1D <- cells_meta %>%
  ggplot(aes(x = factor(sample), group = cell_or_background, fill = cell_or_background, label = after_stat("count")))+
  geom_bar()+
  labs(x = "Sample", y = "SPC Count", fill = NULL)+
  scale_y_continuous(breaks = c(0:10)*100)+
  scale_fill_manual(values = c(cell = "darkblue", background = "grey"))

##plot figure 1E####
figure_1E <- cells_meta %>%
  mutate(sample = paste("Sample", sample)) %>%
  filter(cell_or_background == "cell") %>%
  ggplot(aes(x = n_reads, y = fraction_HU3, color = strain))+
  geom_point()+
  scale_x_continuous(labels = scientific)+
  scale_color_manual(values = c(HU3 = "black", BPK081 = "darkblue", doublet = "darkred", low_cov = "grey"))+
  scale_y_continuous(limits = c(0,1))+
  labs(x = "Total reads", y = "HU3 SNPs/alternative alleles", color = "Strain")+
  facet_wrap(vars(sample), scales = "free_x", nrow = 1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

##plot figure 1F####
figure_1F <- bins_meta %>%
  select(sample, gc_content, mean_raw_counts, mean_corrected_counts) %>%
  filter(sample != "4") %>%
  pivot_longer(cols = c(mean_raw_counts, mean_corrected_counts), values_to = "count", names_to = "count_type") %>%
  mutate(count_type = c(mean_raw_counts = "Raw Counts", mean_corrected_counts = "Corrected Counts")[count_type]) %>%
  group_by(sample, count_type) %>%
  filter(!count %in% boxplot.stats(count)$out) %>%
  mutate(sample = paste("Sample", sample)) %>%
  ggplot(aes(x = gc_content, y = count))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE, color = "blue", linewidth = 1) +  # Linear model line
  scale_x_continuous(limits = c(0.5, 0.7))+
  labs(x = "GC content", y = "Normalized Counts")+
  facet_grid(rows = vars(sample), cols = vars(fct_rev(count_type)), scales = "free")

##plot figure 1G####
#compare raw vs corrected counts
figure_1G <- bins_meta %>%
  group_by(sample) %>%
  mutate(bin_position = row_number()) %>%
  filter(!is_outlier & !is_empty) %>%
  filter(sample != "4") %>%
  mutate(sample = paste("Sample", sample)) %>%
  pivot_longer(cols = c("mean_raw_counts", "mean_corrected_counts"), names_to = "count_type", values_to = "count") %>%
  group_by(sample, chromosome) %>%
  filter(!count %in% boxplot.stats(count)) %>%
  mutate(count_type = c(mean_raw_counts = "Raw Counts", mean_corrected_counts = "Corrected Counts")[count_type]) %>%
  ggplot(aes(x = bin_position, y = count, color = chromosome))+
  geom_point()+
  #geom_boxplot(outliers = FALSE)+
  scale_color_manual(values = rep(c("black", "orange"), times = 100))+
  guides(color = "none")+
  labs(x = "20 kb bin", y = "Mean count")+
  facet_grid(cols = vars(fct_rev(count_type)), rows = vars(sample), scale = "free")

##plot figure 1H#####
##compare lorenz curves with and without GC correction
lorenz_df <- list(raw_counts = lapply(true_cells, function(x)x$counts$raw_counts),
                  corrected_counts = lapply(true_cells, function(x)x$counts$corrected_counts))

lorenz_df <- lapply(lorenz_df, function(x){
  x <- x %>%
    do.call(cbind, .) %>%
    as.data.frame() %>%
    rownames_to_column("bin") %>%
    pivot_longer(cols = -bin, names_to = "cell", values_to = "count") %>%
    mutate(cell = gsub("_.*", "", cell)) %>%
    mutate(sample = true_cells_meta[cell,]$sample) %>%
    group_by(cell) %>%
    arrange(count, .by_group = TRUE) %>%
    mutate(order = row_number()) %>%
    mutate(frac_genome = order/max(order)) %>%
    mutate(frac_count = cumsum(count)) %>%
    mutate(frac_count = frac_count/max(frac_count)) %>%
    group_by(sample, frac_genome) %>%
    summarise(mean = mean(frac_count), se = sd(frac_count)/sqrt(n()))
    return(x)
})

lorenz_df$raw_counts$condition <- "raw"
lorenz_df$corrected_counts$condition <- "corrected"
lorenz_df <- bind_rows(lorenz_df)

lorenz_df %>%
  mutate(condition = fct_rev(factor(condition))) %>%
  ggplot(aes(x = frac_genome, y = mean, color = sample))+
  geom_line()+
  geom_line(data = data.frame(frac_genome = c(0,1), mean = c(0,1)), color = "lightgrey", linetype = "dashed")+
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se), alpha = 0.2, color = NA)+
  facet_wrap(vars(condition))+
  scale_color_manual(values = sample_colors)+
  theme_bw()

##plot figure 1F####
figure_1F <- true_cells_meta %>%
  select(strain, sample, fraction_1, fraction_5, fraction_1_sub) %>%
  pivot_longer(cols = -c(strain, sample), names_to = "fraction") %>%
  mutate(fraction = c(fraction_1 = "covered at least 1x", 
                      fraction_5 = "covered at least 5x", 
                      fraction_1_sub = "covered at least 1x\n(subsampled to 100.000 reads per cell)")[fraction]) %>%
  mutate(fraction = factor(fraction, levels = unique(fraction)[c(1,3,2)])) %>%
  ggplot(aes(x = sample, y = value, fill = sample))+
  geom_violin(scale = "width")+
  geom_boxplot(fill = "white", width = 0.25)+
  scale_fill_manual(values = sample_colors)+
  scale_x_discrete(name = NULL, breaks = NULL)+
  scale_y_continuous(name = "Percentage of the genome")+
  facet_wrap(vars(fraction), nrow = 1, scales = "free")+
  coord_cartesian(clip = FALSE)+
  theme_minimal()

##plot figure 1G####
plot_list <- list()
for(metric in c("gini", "mapd", "ICCV", "ICF_score", "mean_coverage")){
  plot_list[[metric]] <- true_cells_meta %>%
    mutate(plot_label = metric) %>%
    ggplot(aes(x = factor(sample), y = .data[[metric]]))+
    geom_violin(aes(fill = factor(sample)))+
    geom_boxplot(width = 0.075)+
    labs(x = NULL, y = metric)+
    guides(fill = "none")+
    scale_fill_manual(values = sample_colors)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
    facet_wrap(vars(plot_label))
}

#combine them
figure_1G <- patchwork::wrap_plots(plot_list)

#ccombine the figures in a single panel####
top <- figure_1A + figure_1B + plot_layout(widths = c(0.875, 0.125))
middle <- figure_1C + figure_1D + plot_layout(widths = c(0.875, 0.125))
bottom <- figure_1F

final <- (top / middle / figure_1E / figure_1F)+
  plot_annotation(tag_levels = 'A')+ 
  plot_layout(heights = c(0.25, 0.25, 0.25, 0.25))&
  theme(plot.tag = element_text(size = 18, face = "bold"),
    plot.tag.position   = c(0, 1),    # top-left corner
    plot.tag.background = element_rect(fill = NA, colour = NA),
    axis.title = element_text(size = 10),
    text = element_text(size = 10))

#save the final figure panel####
plot_scale <- 1.8
ggsave("figure_1.pdf", plot = final, width = 8.27 * plot_scale, height = 6 * plot_scale)

#open it
if (Sys.info()["sysname"] == "Darwin") {
  system2("open", args = "figure_1.pdf", wait = FALSE)
} else if (Sys.info()["sysname"] == "Linux") {
  system2("xdg-open", args = "figure_1.pdf", wait = FALSE)
}