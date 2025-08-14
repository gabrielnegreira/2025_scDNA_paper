#clean the environment####
source("clean_environment.R")

##to use the new `layout_tags()` function of `ggalign` I had to install it in a different environment... 
##...as its dependency on a beta version of `ggplot2` was breaking `ggtree` in other scripts. Hence, we use...
##... `callr` to call ggalign in that environment instead.
callr::r(function() {
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
  
  #in figure 1 we split the 10X data by strain (as each strain was a different library)
  cells_meta <- cells_meta %>%
    mutate(sample = ifelse(experiment == "10X", paste(sample, strain), sample))
  
  #plot the figures####
  #create an empty list to store the figures
  figures <- list()
  ##plot placeholder for figure 1A####
  
  #Create an empty plot with text and border
  figures[["A"]] <- ggplot(mapping = aes(x = c(0:10), y = c(0:10)))+
    geom_text(aes(label = "Experimental Setup\ndiagram", x = 5, y = 5))+
    theme_void()+
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
  
  ##plot figure 1B####
  figures[["B"]] <- all_SPCs_meta %>%
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
  figures[["C"]] <- all_SPCs_meta %>%
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
  figures[["D"]] <- all_SPCs_meta  %>%
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
  figures[["E"]] <- all_SPCs_meta  %>%
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
  figures[["F"]] <- cells_meta %>%
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
  
    top <- ggalign::align_plots(!!!figures[c("A", "B")], widths = c(0.875, 0.125)) 
    middle <- ggalign::align_plots(!!!figures[c("C", "D")], widths = c(0.875, 0.125))
    final <- ggalign::align_plots(top, middle, !!!figures[c("E", "F")], ncol = 1)
    
    #add layout tags
    final <- final + layout_tags("A") & theme(plot.tag = element_text(size = 16))
    
    #save the final figure panel####
    plot_scale <- 1.8
    ggsave("figure_1.pdf", plot = final, width = 8.27 * plot_scale, height = 7 * plot_scale)
}, show = TRUE)

#open it
if (Sys.info()["sysname"] == "Darwin") {
  system2("open", args = "figure_1.pdf", wait = FALSE)
} else if (Sys.info()["sysname"] == "Linux") {
  system2("xdg-open", args = "figure_1.pdf", wait = FALSE)
}

