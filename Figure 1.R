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
all_SPCs <- readRDS("inputs/karyotyping_objects/all_SPCs.rds")
cells_meta <- read_delim("inputs/cell_qc/cells_meta.tsv") %>%
  column_to_rownames("rowname")

all_SPCs_meta <- lapply(all_SPCs, function(x)x$metadata$cells_meta) %>%
  bind_rows()
##set base theme parameters for ggplot####
update_geom_defaults("point", list(size = 0.5))

#in figure 1 we split the 10X data by strain (as each strain was a different library)
cells_meta <- cells_meta %>%
  mutate(sample = ifelse(experiment == "10X", paste(sample, strain), sample))

#plot the figures####
##plot placeholder for figure 1A####

#Create an empty plot with text and border
figure_1A <- ggplot()+
  annotate("text", x = 0.5, y = 0.5, label = "Experimental Setup\nSchematic", size = 6, hjust = 0.5)+
  theme_void()+
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE)+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

##plot figure 1B####
figure_1B <- all_SPCs_meta %>%
  filter(experiment == "Atrandi") %>%
  ggplot(aes(x = n_reads, y = fct_rev(factor(sample)), fill = factor(sample)))+
  geom_density_ridges(alpha = 0.5)+
  scale_x_continuous(transform = "log10", breaks = c(2000, 17000, 170000, 900000))+
  scale_y_discrete(name = NULL)+
  scale_fill_manual(values = sample_colors)+
  labs(x = "Total reads per cell (log10 scale)", y = "Sample")+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid.minor = element_blank())

##plot figure 1C####
figure_1C <- all_SPCs_meta %>%
  filter(experiment == "Atrandi") %>%
  group_by(sample) %>%
  arrange(desc(n_reads)) %>%
  mutate(rank = row_number()) %>%
  group_by(sample, cell_or_background) %>%
  mutate(last_true_cell = rank == rank[which(rank == max(rank))]) %>%
  ungroup() %>%
  mutate(last_true_cell = last_true_cell & cell_or_background == "cell") %>%
  mutate(label = ifelse(last_true_cell, paste0("Rank: ", rank, "\nread_count: ", n_reads), NA)) %>%
  ggplot(aes(x = rank, y = n_reads, color = cell_or_background, label = label))+
  geom_line()+
  geom_text_repel(show.legend = FALSE, nudge_x = -1, nudge_y = -0.5, size = 2, min.segment.length = 0)+
  scale_x_continuous(transform = "log10", limits = c(1, NA), breaks = c(1, 10, 100, 1000))+
  scale_y_continuous(transform = "log10")+
  scale_color_manual(values = c(cell = "darkblue", backgroun = "lightgrey"))+
  guides(color = "none")+
  labs(x = "Barcode Rank (log10)", y = "Read Count (log10)")+
  facet_wrap(vars(sample), nrow = 1, scales = "free")+
  theme_minimal()+
  theme(panel.grid.minor = element_blank())

##plot figure 1D####
#check how many cells were removed from each sample
figure_1D <- all_SPCs_meta  %>%
  filter(experiment == "Atrandi") %>%
  ggplot(aes(x = factor(sample), group = cell_or_background, fill = cell_or_background, label = after_stat("count")))+
  geom_bar()+
  labs(x = NULL, y = "SPC Count", fill = NULL)+
  scale_y_continuous(breaks = c(0:10)*100)+
  scale_fill_manual(values = c(cell = "darkblue", background = "grey"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#remove background cells####
all_SPCs_meta <- filter(all_SPCs_meta, cell_or_background == "cell" & sample != "Sample 4")

##plot figure 1E####
figure_1E <- all_SPCs_meta  %>%
  filter(experiment == "Atrandi") %>%
  ggplot(aes(x = n_reads, y = fraction_HU3, color = strain))+
  geom_point()+
  scale_x_continuous(labels = scientific)+
  scale_color_manual(values = c(HU3 = "black", BPK081 = "darkblue", doublet = "darkred", low_cov = "grey"))+
  scale_y_continuous(limits = c(0,1))+
  labs(x = "Total reads", y = "HU3 SNPs/alternative alleles", color = "Strain")+
  facet_wrap(vars(sample), scales = "free_x", nrow = 1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

##remove doublets####
cells_meta <- filter(cells_meta, strain != "doublet")

##plot figure 1F####
figure_1F <- cells_meta %>%
  #filter(experiment == "Atrandi") %>%
  select(strain, sample, fraction_1, fraction_5, fraction_1_sub, mean_coverage) %>%
  pivot_longer(cols = -c(strain, sample), names_to = "fraction") %>%
  mutate(fraction = c(mean_coverage = "Mean coverage (X)", 
                      fraction_1 = "fraction of the genome covered at least 1x", 
                      fraction_5 = "fraction of the genome covered at least 5x", 
                      fraction_1_sub = "fraction of the genome covered at least 1x\n(subsampled to 100.000 reads per cell)")[fraction]) %>%
  mutate(fraction = factor(fraction, levels = unique(fraction)[c(4,1,3,2)])) %>%
  ggplot(aes(x = sample, y = value, fill = sample))+
  geom_violin(scale = "width")+
  geom_boxplot(fill = "white", width = 0.25, outlier.size = 0.5)+
  scale_fill_manual(values = sample_colors)+
  scale_x_discrete(name = NULL)+
  scale_y_continuous(name = NULL)+
  facet_wrap(vars(fraction), nrow = 1, scales = "free")+
  guides(fill = NULL)+
  coord_cartesian(clip = FALSE)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())


#combine the figures in a single panel####
top <- figure_1A + figure_1B + plot_layout(widths = c(0.875, 0.125))
middle <- figure_1C + figure_1D + plot_layout(widths = c(0.875, 0.125))

final <- (top / middle / figure_1E / figure_1F)+
  plot_annotation(tag_levels = 'A')+ 
  plot_layout(heights = c(0.25, 0.20, 0.3, 0.3))&
  theme(plot.tag = element_text(size = 18, face = "bold"),
    plot.tag.position   = c(0, 1),    # top-left corner
    plot.tag.background = element_rect(fill = NA, colour = NA),
    axis.title = element_text(size = 10),
    text = element_text(size = 10))

#save the final figure panel####
plot_scale <- 1.8
ggsave("figure_1.pdf", plot = final, width = 8.27 * plot_scale, height = 7 * plot_scale)

#open it
if (Sys.info()["sysname"] == "Darwin") {
  system2("open", args = "figure_1.pdf", wait = FALSE)
} else if (Sys.info()["sysname"] == "Linux") {
  system2("xdg-open", args = "figure_1.pdf", wait = FALSE)
}