#clean the environment####
rm(list = ls(envir = parent.frame()))
gc()

#load libraries####
library(tidyverse)
library(ggridges)
library(ggrepel)
library(Cairo)
library(cowplot)
library(patchwork)
library(scales)

#get inputs####
scDNAobj_list <- readRDS("inputs/true_cells.rds")

##set base theme parameters for ggplot####
update_geom_defaults("point", list(size = 0.5))

#bind cells meta####
cells_meta <- list()
for(i in c(1:length(scDNAobj_list))){
  cells_meta[[i]] <- scDNAobj_list[[i]]$metadata$cells_meta
}
cells_meta <- bind_rows(cells_meta)

#bind bins_meta####
bins_meta <- list()
for(i in c(1:length(scDNAobj_list))){
  bins_meta[[i]] <- scDNAobj_list[[i]]$metadata$bins_meta
  bins_meta[[i]]$sample <- unique(scDNAobj_list[[i]]$metadata$cells_meta$sample)
}
bins_meta <- bind_rows(bins_meta)

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
  scale_x_continuous(transform = "log10")+
  labs(x = "Total read count", y = "Sample")+
  guides(fill = "none")

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
  facet_grid(rows = vars(sample), cols = vars(fct_rev(count_type)), scales = "free")#+
  #theme(axis.text = element_text(size = 14), axis.title = element_text(size = 18), strip.text = element_text(size = 14))

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
  facet_grid(cols = vars(fct_rev(count_type)), rows = vars(sample), scale = "free")#+
  #theme(axis.text = element_text(size = 14), axis.title = element_text(size = 18), strip.text = element_text(size = 14))

#ccombine the figures in a single panel####
top <- figure_1A + figure_1B + plot_layout(widths = c(0.875, 0.125))
middle <- figure_1C + figure_1D + plot_layout(widths = c(0.875, 0.125))
bottom <- figure_1F + figure_1G + plot_layout(widths = c(0.3, 0.7))

final <- (top / middle / figure_1E / bottom)+
  plot_annotation(tag_levels = 'A')+ 
  plot_layout(heights = c(0.125, 0.125, 0.15, 0.7))&
  theme(plot.tag = element_text(size = 18, face = "bold"),
    plot.tag.position   = c(0, 1),    # top-left corner
    plot.tag.background = element_rect(fill = NA, colour = NA),
    axis.title = element_text(size = 10),
    text = element_text(size = 10))

#save the final figure
plot_scale <- 1.8
ggsave("figure_1.pdf", plot = final, width = 8.27 * plot_scale, height = 9 * plot_scale)

#open it
system2('open', args = "figure_1.pdf", wait = FALSE)