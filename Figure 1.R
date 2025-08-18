#clean the environment####
source("clean_environment.R")

##to use the new `layout_tags()` function of `ggalign` I had to install it in a different environment... 
##...as its dependency on a beta version of `ggplot2` was breaking `ggtree` in other scripts. Hence, we use...
##... `libPaths` to call ggalign in that environment instead.

.libPaths("~/R/lib-ggalign_beta")
#load libraries####
library(tidyverse)
library(ggridges)
library(ggrepel)
library(ggalign)
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

#plot the figures####
#create an empty list to store the figures
figures <- list()

##plot placeholder for figure 1A####
#Create an empty plot with text and border
figures[[length(figures) + 1]] <- ggplot(mapping = aes(x = c(0:10), y = c(0:10)))+
  geom_text(aes(label = "Experimental Setup\ndiagram", x = 5, y = 5))+
  theme_void()+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

##plot figure 1B####
figures[[length(figures) + 1]] <- all_SPCs_meta %>%
  filter(experiment == "Atrandi") %>%
  ggplot(aes(x = n_reads, y = fct_rev(factor(sample)), fill = factor(sample)))+
  geom_density_ridges(alpha = 0.5)+
  scale_x_continuous(transform = "log10", breaks = c(2000, 17000, 170000, 900000))+
  scale_y_discrete(name = NULL)+
  scale_fill_manual(values = sample_colors)+
  labs(x = "Total reads per cell (log10 scale)", y = "Sample")+
  guides(fill = "none")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid.minor = element_blank())

##plot figure 1C####
figures[[length(figures) + 1]] <- all_SPCs_meta %>%
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
  theme_bw()+
  theme(panel.grid.minor = element_blank())

##plot figure 1D####
#check how many cells were removed from each sample
figures[[length(figures) + 1]] <- all_SPCs_meta  %>%
  filter(experiment == "Atrandi") %>%
  ggplot(aes(x = factor(sample), group = cell_or_background, fill = cell_or_background, label = after_stat("count")))+
  geom_bar()+
  labs(x = NULL, y = "SPC Count", fill = NULL)+
  scale_y_continuous(breaks = c(0:10)*100)+
  scale_fill_manual(values = c(cell = "darkblue", background = "grey"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid.minor = element_blank())

#remove background cells####
all_SPCs_meta <- filter(all_SPCs_meta, cell_or_background == "cell" & sample != "Sample 4")

##plot figure 1E####
figures[[length(figures) + 1]] <- all_SPCs_meta  %>%
  filter(experiment == "Atrandi") %>%
  ggplot(aes(x = n_reads, y = fraction_HU3, color = strain))+
  geom_point()+
  scale_x_continuous(labels = scientific)+
  scale_color_manual(values = c(HU3 = "black", BPK081 = "darkblue", doublet = "darkred", low_cov = "grey"))+
  scale_y_continuous(limits = c(0,1))+
  labs(x = "Total reads", y = "fraction HU3 signature", color = "Strain")+
  facet_wrap(vars(sample), scales = "free_x", nrow = 1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.minor = element_blank())

##remove doublets####
cells_meta <- filter(cells_meta, strain != "doublet")

##plot figure 1F####
for(fraction in c("mean_coverage", "fraction_1", "fraction_5", "fraction_1_sub")){
  figures[[length(figures) + 1]] <- cells_meta %>%
    ggplot(aes(x = sample, y = .data[[fraction]], fill = sample))+
    geom_violin(scale = "width")+
    geom_boxplot(fill = "white", width = 0.25, outlier.size = 0.5)+
    scale_fill_manual(values = sample_colors)+
    scale_x_discrete(name = NULL)+
    scale_y_continuous(name = c(`mean_coverage` = "coverage (X)",
                                `fraction_1` = "genome fraction",
                                `fraction_5` = "genome fraction",
                                `fraction_1_sub` = "genome fraction")[fraction])+
    guides(fill = "none")+
    ggtitle(c(`mean_coverage` = "Mean Coverage (X)",
              `fraction_1` = "Fraction of genome covered at least 1X",
              `fraction_5` = "Fraction of genome covered at least 5X",
              `fraction_1_sub` = "Fraction of genome covered at least 1X\n(Subsampled to 100.000 reads per cell)")[fraction])+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 8))
}

  #combine the figures in a single panel####
  top <- ggalign::align_plots(!!!figures[c(1,2)], widths = c(0.875, 0.125)) 
  middle <- ggalign::align_plots(free_border(figures[[3]], borders = "b"), figures[[4]], widths = c(0.875, 0.125), )
  bottom <- ggalign::align_plots(!!!figures[c(6:length(figures))], nrow = 1)
  final <- ggalign::align_plots(top, middle, figures[[5]], bottom, 
                                ncol = 1, 
                                heights = c(0.3, 0.2, 0.25, 0.25))
  
  #add layout tags
  final <- final + layout_tags("A") + layout_theme(plot.tag = element_text(size = 16))
  
  #save the final figure panel####
  plot_scale <- 1.8
  ggsave("figure_1.pdf", plot = final, width = 8.27 * plot_scale, height = 6.5 * plot_scale)

#open it
if (Sys.info()["sysname"] == "Darwin") {
  system2("open", args = "figure_1.pdf", wait = FALSE)
} else if (Sys.info()["sysname"] == "Linux") {
  system2("xdg-open", args = "figure_1.pdf", wait = FALSE)
}

