#clean the environment####
rm(list = ls(envir = parent.frame()))
gc()

#load libraries####
source("color_palettes.R")
library(tidyverse)
library(ggridges)
library(ggalign)

#get inputs
cells_meta <- read_delim("inputs/cell_qc/cells_meta.tsv") %>%
  column_to_rownames("rowname") %>%
  filter(cell_or_background == "cell" & strain != "doublet" & sample != "Sample 4")
samples <- c(`sample_1` = "Sample 1", `sample_2` = "Sample 2", `sample_3` = "Sample 3", `sample_5` = "Sample 5", `sample_6` = "Sample 6")
samples_to_plot <- c("Sample 2", "Sample 6")
chromosomes <- c("Ld29", "Ld36")
regions <- list(Ld29 = c(9613, 9647), Ld36 =c(25620, 25679))

#create a key/pair vector containing the correct barcode string for each barcode
bcs <- cells_meta %>%
  rownames_to_column("full_barcode") %>%
  mutate(partial_barcode = gsub("_.*", "", full_barcode)) %>%
  select(partial_barcode, full_barcode) %>%
  deframe()

#create a key/value pair for cells strains
strains <- cells_meta %>%
  rownames_to_column("cell") %>%
  select(cell, strain) %>%
  deframe()


#get the cnv files
cnv_list <- list()
hm_list <- list()
dens_list <- list()
cells_meta[chromosomes] <- NA
for(Sample in names(samples)){
  cnv_list[[Sample]] <- list()
  sample_name <- samples[Sample]
  for(Chromosome in chromosomes) {
    Pattern <- paste0("^", Sample, ".*", Chromosome, "_100nt_CNV\\.csv$")
    cnv_file <-  list.files("inputs/cnvs/", pattern = Pattern)
    cnv_file <- paste0("inputs/cnvs/", cnv_file)

    mat <- read.csv(cnv_file, row.names = 1) %>%
      as.data.frame() %>%
      mutate(rowname = paste0(chrom, "_", start)) %>%
      column_to_rownames("rowname") %>%
      select(-chrom, -start, -end) %>%
      as.matrix()
    
    #fix matrix colnames
    colnames(mat) <- gsub("_", "", colnames(mat))
    colnames(mat) <- bcs[colnames(mat)]
    
    #remove doublets
    cells <- names(strains)
    cells <- cells[which(cells %in% colnames(mat))]
    mat <- mat[,cells]
    
    #covert infinite values to 0
    mat[is.infinite(mat)] <- 0

    #normalize the matrix
    mat[,] <- apply(mat, 2, normalize, method = "mean")
    
    #save the matrix
    cnv_list[[Sample]][[Chromosome]] <- mat
    
    #estimate copy number for the cnv and save it in the cells_meta
    indices <- regions[[Chromosome]]
    indices <- indices[1]:indices[2]
    copy_numbers <- colMeans(mat[indices,])/colMeans(mat[-indices,])
    cells_meta <- cells_meta %>%
      rownames_to_column("rowname") %>%
      mutate(!!Chromosome := ifelse(is.na(.data[[Chromosome]]), copy_numbers[rowname], .data[[Chromosome]])) %>%
      column_to_rownames("rowname")
      
    
    if(sample_name %in% samples_to_plot){
      #trim the matrix to focus on the locus
      indices <- regions[[Chromosome]]
      indices[1] <- indices[1] - 100
      indices[2] <- indices[2] + 100
      indices <- indices[1]:indices[2]
      mat <- mat[indices,]
      rownames(mat) <- indices*100
      #set row breaks
      breaks <- round(seq(1, nrow(mat), by = 10))
      breaks <- rownames(mat)[breaks]
      breaks <- setNames(as.integer(breaks), breaks)
      breaks <- format(breaks, big.mark = ",", scientific = FALSE)
      #set color scale
      color_range <- quantile(mat, c(0, 0.99))
      hm <- ggheatmap(mat)+
        scale_x_discrete(breaks = NULL)+
        scale_y_continuous(breaks = names(breaks), labels = breaks)+
        theme(axis.text = element_text())+
        scale_fill_scico(palette = "lipari", limits = color_range, na.value = scico(100, palette = "lipari")[100])+
        anno_top(size = 0.3)+
        ggalign(size = 0)+
        patch_titles(top = paste0(sample_name, "; Chromosome: ", Chromosome))+
        align_dendro()+
        theme_void()+
        anno_top() + #create a second annotation space, below the dendogram and on top of the heatmap
        ggalign(data = strains[colnames(mat)], size = 0.2)+ #initiate a ggplot object withy the strains as data.
        geom_tile(aes(y = 1, fill = factor(value))) + #add annotation tiles.
        scale_fill_manual(name = "Strain", values = strain_colors)+ #set the annotation colors
        theme_void()
      
      
      hm_list[[paste0(Sample, ": ", Chromosome)]] <- hm
      
    }
  }
}


#plot density
density_plot <- cells_meta %>%
  filter(experiment == "Atrandi") %>%
  select(sample, strain, matches(chromosomes)) %>%
  pivot_longer(cols = -c(sample, strain), names_to = "Chromosome", values_to = "copy_number") %>%
  mutate(CNV = paste("CNV in chromosome", Chromosome)) %>%
  ggplot(aes(x = copy_number, fill = strain, group = strain))+
  geom_density(alpha = 0.5)+
  scale_fill_manual(values = strain_colors)+
  scale_x_continuous(n.breaks = 10) +
  labs(x = "haploid copy number", fill = "strain", y = "density")+
  facet_grid(rows = vars(CNV), cols = vars(sample), scales = "free")+
  theme_bw()+
  theme(panel.grid.minor = element_blank())

hm_plot <- ggalign::align_plots(!!!hm_list)

final <- align_plots(hm_plot, NULL, density_plot, ncol = 1, heights = c(0.69, 0.02, 0.29))

#save the final figure panel####
plot_scale <- 1.8
ggsave("figure_4.pdf", plot = final, width = 8.27 * plot_scale, height = 9 * plot_scale)

#open it
if (Sys.info()["sysname"] == "Darwin") {
  system2("open", args = "figure_4.pdf", wait = FALSE)
} else if (Sys.info()["sysname"] == "Linux") {
  system2("xdg-open", args = "figure_4.pdf", wait = FALSE)
}
