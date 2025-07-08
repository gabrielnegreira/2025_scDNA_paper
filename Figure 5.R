#clean the environment####
rm(list = ls(envir = parent.frame()))
gc()

#load libraries####
library(tidyverse)
library(patchwork)
library(vcfR)
library(xlsx)
library(fuzzyjoin)
library(ggalign)
source("color_palettes.R")


#get inputs####
## get scDNA object of true cells
true_cells <- read_rds("inputs/true_cells.rds")

#plot figure 5A####
##get PCA eigen vectors
pca_files <- list.files("inputs", pattern = "pca.eigenvec")
pca_data <- list()
for(i in seq_along(pca_files)){
  Sample <-  as.integer(gsub(".*sample_(\\d+)_.*", "\\1", pca_files[[i]]))
  index = ifelse(Sample >= 5, Sample - 1, Sample)
  pca_data[[i]] <- read_delim(paste0("inputs/", pca_files[[i]])) %>%
    select(-`#FID`) %>%
    rename(cell = IID) %>%
    mutate(cell = gsub("_", "", cell)) %>%
    mutate(sample = Sample) %>%
    mutate(strain = true_cells[[index]]$metadata$cells_meta[cell,]$strain)
}

##merge them in a single data frame
pca_data <- bind_rows(pca_data)
##plot it
figure_5A <- pca_data %>%
  mutate(sample = paste("Sample", sample)) %>%
  mutate(strain = ifelse(is.na(strain), "doublet", strain)) %>%
  ggplot(aes(x = PC1, y = PC2, color = strain))+
  geom_point(size = 0.5)+
  facet_wrap(vars(sample), scales = 'free', nrow = 1)+
  guides(color = guide_legend(override.aes = list(size = 2)))+
  scale_color_manual(values = c(strain_colors, doublet = "grey"))+
  theme_bw()


#plot figure 5B####
##get drug resistance snps
drug_data <- read.xlsx("inputs/drug_resistance_Ldon_v2.xlsx", sheetIndex = 1)
plot_list <- list()
for(Sample in c(5, 6)){
  #get needed files
  vcf <-  read.vcfR(paste0("inputs/sample_", Sample, ".filtered.drugs.snpeff.vcf"))
  snpeff_data <- read.delim(paste0("inputs/sample_", Sample, ".filtered.drugs.snpeff.csv"))
  #add the drug_data to the sneff_data
  snpeff_data <- fuzzy_inner_join(snpeff_data, drug_data,
                                  by = c("CHROM" = "chrom", "POS" = "start", "POS" = "end"),
                                  match_fun = list(`==`, `>=`, `<=`))
  
  #remove non impactful SNPs
  info_data <- INFO2df(vcf)
  ann_data <- str_split_fixed(info_data$ANN, "\\|", n=10)
  variants_omit <- c("start_lost", "downstream_gene_variant", "synonymous_variant", "upstream_gene_variant", "intergenic_region")
  vcf <- vcf[! ann_data[,2] %in% variants_omit, ]
  
  #convert the vcf file to a matrix
  x <- vcfR2genlight(vcf)
  gt_matrix = t(as.data.frame(x))
  
  ## remove those SNPs that are the same (or NA) in all cells
  gt_matrix <- gt_matrix[apply(gt_matrix, 1, function(row) {
    unique_values <- unique(na.omit(row)) # Remove NA values and get unique values
    length(unique_values) > 1            # Keep rows with more than one unique value
  }), ]
  
  ##get cells' strain
  ###match cell barcodes to the metadata
  colnames(gt_matrix) <- gsub("_", "", colnames(gt_matrix))
  index <- ifelse(Sample >= 5, Sample - 1, Sample)
  cells_strain <- true_cells[[index]]$metadata$cells_meta %>%
    rownames_to_column("cell") %>%
    mutate(cell = gsub("_.*", "", cell)) %>%
    dplyr::select(cell, strain) %>%
    deframe()
  cells_strain <- cells_strain[colnames(gt_matrix)]
  names(cells_strain) <- colnames(gt_matrix) #cells with NA values lose their name. Fix it.
  cells_strain[is.na(cells_strain)] <- "doublet"
  
  #remove doublets
  gt_matrix <- gt_matrix[,names(cells_strain[cells_strain != "doublet"])]
  
  ## group genes by drug
  gene_groups <- snpeff_data %>%
    mutate(row = paste0(CHROM, "_", POS)) %>%
    mutate(Drug_code = c(MSL = "MIL", MIL = "MIL", Hlocus = "ANT", AQP = "ANT", Amb = "AMB")[Drug_code]) %>%
    mutate(Gene_code = gsub("AQP", "AQP1", Gene_code)) %>%
    mutate(label = paste0(Drug_code, ": ", Gene_code)) %>%
    dplyr::select(row, label) %>%
    deframe()
  
  #update row names
  gene_labels <- snpeff_data %>%
    mutate(row = paste0(CHROM, "_", POS)) %>%
    mutate(label = paste0(row, "::", ANN.0..HGVS_P)) %>%
    dplyr::select(row, label) %>%
    deframe()
  
  rownames(gt_matrix) <- gene_labels[rownames(gt_matrix)]
  names(gene_groups) <- gene_labels[names(gene_groups)]
  
  #convert the gt_matrix to integer values
  
  #plot the heatmap
  #to do: swap row and row group label sides
  #       remove AQP1 frame shift mutation in sample 6
  #       highlight SNPs common in both datasets
  
  gt_matrix[,] <- apply(gt_matrix, 2, factor)
  plot <- ggheatmap(gt_matrix)+
    scale_fill_manual(name = "Allele", breaks = c(0, 1, 2), labels = c(`0` = "Homozygous reference", `1` = "Heterozygous", `2` = "Homozygous alternative"), values = snp_colors)+
    scale_x_discrete(breaks = NULL, name = NULL)+
    #facet_grid(switch = "y")+ #this will place the row group labels on the left side instead of the right side. 
    theme(axis.text.y = element_text(size = 6),
          strip.text.y.right = element_text(angle = 0, hjust = 0, vjust = 0.5, size = 6), #this is needed to keep panel labels (gene groups) as default `theme_ggalign()` removes them. 
          strip.placement = "outside", #allow labels to be draw outside the box
          strip.clip = "off") + #prevent labels from being clipped. 
    anno_top(size = 0.1)+
    align_dendro(size = 0.1)+
    patch_titles(top = paste("Sample", Sample))+ #add a title to the top annotation
    theme_void()+
    anno_left()+
    align_group(group = gene_groups[rownames(gt_matrix)])+
    anno_top()+
    ggalign(data = cells_strain[colnames(gt_matrix)], size = 0.1)+ #initiate a ggplot object withy the strains as data.
    geom_tile(aes(y = 1, fill = factor(value))) + #add annotation tiles.
    scale_fill_manual(name = "Strain", values = c(strain_colors)) + #set the annotation colors
    scale_y_discrete(name = "Strain")+
    theme_void()+
    theme(axis.title.y = element_text(size = 8))&
    theme(
      #legend.position = "bottom",  # << move to bottom
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 7),
      legend.key.size = unit(0.3, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0)
    )
 
  plot_list <- c(plot_list, plot)
  
}

figure_5B <- ggalign::align_plots(!!!plot_list, guides = "rlb")
figure_5B <- patchwork::wrap_plots(ggalignGrob(figure_5B), ncol = 1) # compactly stack

final <- figure_5A + figure_5B + 
  plot_layout(heights = c(0.2, 0.8))+
  plot_annotation(tag_levels = 'A')&
  theme(plot.tag = element_text(size = 18, face = "bold"),
        plot.tag.position   = c(0, 1),    # top-left corner
        plot.tag.background = element_rect(fill = NA, colour = NA),
        axis.title = element_text(size = 10),
        text = element_text(size = 10),
        plot.margin = margin(5, 5, 5, 5))

#save the final figure panel####
plot_scale <- 1.8
ggsave("figure_5.pdf", plot = final, width = 8.27 * plot_scale, height = 6 * plot_scale)

#open it
if (Sys.info()["sysname"] == "Darwin") {
  system2("open", args = "figure_5.pdf", wait = FALSE)
} else if (Sys.info()["sysname"] == "Linux") {
  system2("xdg-open", args = "figure_5.pdf", wait = FALSE)
}

