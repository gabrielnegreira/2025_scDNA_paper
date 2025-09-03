#libraries
library(tidyverse)
library(ggalign)
library(readr)
library(pracma)
library(mgcv)

#parameters
ladder_sizes <- c(15, 100, 250, 400, 600, 1000, 1500, 2500, 3500, 5000, 10000)

#get inputs
file_names <- list.files("inputs/tapestation/")
files <- lapply(paste0("inputs/tapestation/", file_names), read_csv)
names(files) <- file_names

#prepare files for plotting
files <- lapply(files, function(x){
  ## convert to long dataframe
  x <- x %>%
    rownames_to_column("time") %>%
    mutate(time = as.integer(time)) %>%
    pivot_longer(names_to = "sample", values_to = "fluorescence", cols = -time) %>%
    mutate(sample_levels = as.integer(as.factor(sample))) %>%
    mutate(sample = gsub("*.1: ", "", sample),
           sample = gsub("Sample 1", "10 pg Human DNA", sample),
           sample = gsub("Sample 2", "0.1 pg Human DNA", sample),
           sample = gsub("Sample 3", "10 pg Leishmania DNA", sample),
           sample = gsub("Sample 4", "0.1 pg Leishmania DNA", sample),
           sample = gsub("Sample 5", "Negative Control", sample),
           sample = gsub(" - 1 in 10", "*", sample),
           sample = gsub("Leish ", "Leishmania ", sample),
           sample = gsub("100 pg Leishmania DNA*", "Positive Control", sample),
           sample = gsub("EB Buffer", "Negative Control", sample),
           sample = gsub("#.*- ", "", sample),
           sample = fct_reorder(sample, sample_levels))
})

# find the peaks in the ladded sample
ladder_trace <- files$pta_test1.csv %>% 
  filter(sample == "Ladder")
peaks <- findpeaks(ladder_trace$fluorescence, minpeakheight = 200, minpeakdistance = 10)
peaks <- data.frame(time = sort(peaks[,2]), size = ladder_sizes)

plot_list <- lapply(files, function(x){
  x %>%
    filter(sample != "Ladder") %>%
    ggplot(aes(x = time, y = fluorescence))+
    geom_line()+
    geom_area(fill = "lightgrey", alpha = 0.5)+
    facet_wrap(vars(sample), ncol = 1, scales = "free")+
    theme_bw()+
    scale_x_continuous(name = "size (bp)", breaks = peaks$time, labels = peaks$size)+
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
})

final <- ggalign::align_plots(!!!plot_list, nrow = 1) + layout_tags("A")
final <- final + layout_theme(plot.tag = element_text(size = 16))

#save the final figure panel####
plot_scale <- 1.8
ggsave("supplementary_figure_2.pdf", plot = final, width = 8.27 * plot_scale, height = 6.5 * plot_scale)

#open it
if (Sys.info()["sysname"] == "Darwin") {
  system2("open", args = "supplementary_figure_2.pdf", wait = FALSE)
} else if (Sys.info()["sysname"] == "Linux") {
  system2("xdg-open", args = "supplementary_figure_2.pdf", wait = FALSE)
}