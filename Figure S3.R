#clean the environment####
source("clean_environment.R")

#load libraries####
library(tidyverse)
library(ggalign)
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
update_geom_defaults("point", list(size = 0.001, shape = 1))
update_geom_defaults("boxplot", list(outlier.size = 0.5))

#fix sample names
names(all_SPCs) <- sample_names[names(all_SPCs)]
names(true_cells) <- sample_names[names(true_cells)]

#plot the figures####
#create an empty list to store the figures
figures <- list()


##plot effect of gc on counts####
#prepare the data for plotting
to_plot <- list("raw_counts", "corrected_counts")
to_plot <- lapply(
  to_plot,
  function(type){
    true_cells %>%
      lapply(
        function(x){
          x <- x$counts[[type]] %>%
            as.data.frame() %>%
            rownames_to_column("bin") %>%
            pivot_longer(cols = -bin, names_to = "barcode", values_to = "count") %>%
            left_join(
              x$metadata$cells_meta %>%
                select(barcode, sample, strain),
              by = "barcode"
            ) %>%
            left_join(
              x$metadata$bins_meta %>%
                mutate(bin_number = row_number()) %>%
                select(bin, bin_number, chromosome, gc_content),
              by = "bin"
            ) %>%
            mutate(count_type = type) %>%
            group_by(barcode) %>%
            mutate(count = normalize(count))
        }
      ) %>%
      bind_rows(.id = "sample")
  }
) %>%
bind_rows() %>%
  group_by(sample, strain, count_type, chromosome, bin_number) %>%
  summarise(mean_count = mean(count), gc_content = mean(gc_content)) %>%
  mutate(sample = factor(sample, levels = sample_names)) %>% 
  mutate(count_type = factor(count_type, levels = c("raw_counts", "corrected_counts")))


#plot the graphs
plot_list <- list()
index <- 0
for(Strain in unique(to_plot$strain)){
  index <- index + 1
  plot_list[[index]] <- to_plot %>%
    filter(strain == Strain) %>%
    ggplot(aes(x = gc_content, y = mean_count))+
    geom_point()+
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5)+
    labs(title = Strain, x = "GC content", y = "normalized mean count")+
    facet_grid(rows = vars(count_type), cols = vars(sample))+
    scale_x_continuous(limits = c(0.5, 0.7))+
    scale_y_continuous(limits = c(0, 3))
  
  index <- index + 1
  plot_list[[index]] <- to_plot %>%
    ggplot(aes(x = bin_number, y = mean_count, color = chromosome))+
    geom_point()+
    guides(color = "none")+
    scale_color_manual(values = rep(c("#719bb7", "#004472"), times = 100))+
    labs(title = Strain, x = "20kb bin", y = "normalized mean count")+
    facet_grid(cols = vars(count_type), rows = vars(sample))+
    scale_y_continuous(limits = c(0, 3))
}

top <- ggalign::align_plots(plot_list[[1]], plot_list[[3]]) + layout_tags(NULL)
bottom <- ggalign::align_plots(plot_list[[2]], plot_list[[4]]) + layout_tags(NULL)
final <- ggalign::align_plots(top, bottom, ncol = 1, heights = c(0.3, 0.6)) + layout_tags("A")
final

#save the final figure panel####
plot_scale <- 1.8
ggsave("figure_S3.pdf", plot = final, width = 8.27 * plot_scale, height = 7 * plot_scale)


#open it
if (Sys.info()["sysname"] == "Darwin") {
  system2("open", args = "figure_S3.pdf", wait = FALSE)
} else if (Sys.info()["sysname"] == "Linux") {
  system2("xdg-open", args = "figure_S3.pdf", wait = FALSE)
}

