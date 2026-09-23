#clean the environment####
source("clean_environment.R")

#load libraries####
library(tidyr)
library(dplyr)
library(tibble)
library(readr)
library(ggalign)
library(patchwork)
library(vcfR)
library(cluster)
library(ape)
library(xlsx)
library(fuzzyjoin)
library(ggtree)
library(ggstar)
library(parallel)
source("color_palettes.R")

#set inputs path
eig_vectors_dir <- "inputs/nucleotide_variants/plink_pca/"
drugs_data_dir <- "inputs/nucleotide_variants/drugs_data/"
allele_depth_dir <- "inputs/nucleotide_variants/allele_depth_matrices/"
scDNAobjs_dir <- "inputs/karyotyping_objects/"
cell_meta_dir <- "inputs/cell_qc/"

#get inputs####
## get scDNA objects
true_cells <- read_rds(paste0(scDNAobjs_dir, "/true_cells.rds"))
## get genomic sites of drug resistance markers
drug_data <- read.xlsx(paste0(drugs_data_dir, "/drug_resistance_Ldon_v2.xlsx"), sheetIndex = 1)
#get cells metadata for all SPCs
cells_meta <- read_rds(paste0(scDNAobjs_dir, "/all_SPCs.rds")) %>%
  lapply(function(x) x$metadata$cells_meta) %>%
  bind_rows() %>%
  mutate(
    sample = sample_names[sample],
    simple_barcode = gsub("_.*", "", barcode) #since pieter did the pca data without the numeric index in the barcode, I create an extra column with barcodes without this index.
    )
  
#plot figure 5A####
##get PCA eigen vectors and eigen values
### first adjust sample_names to match the pca files
names(sample_names) <- gsub("Sample ", "sample_", names(sample_names)) #make the sample names match the pca files
sample_names <- sample_names[grepl("sample_", names(sample_names))] #keep only the spc samples
sample_names <- sample_names[names(sample_names) != "sample_4"]

#now find the files
eig_vectors_files <- list.files(eig_vectors_dir, pattern = "pca.eigenvec")
eig_values_files <- list.files(eig_vectors_dir, pattern = "pca.eigenval")

#load the files and format them for each sample
eig_vectors <- list()
eig_values <- list()

for(Sample in names(sample_names)){
  #get the file names for the sample
  eig_vec <- eig_vectors_files[grepl(Sample, eig_vectors_files)]
  eig_val <- eig_values_files[grepl(Sample, eig_values_files)]
  
  #read the files
  ##get the eigen values and store in a list
  eig_vec <-  read_delim(paste0(eig_vectors_dir, eig_vec))
  PCs <- colnames(eig_vec)[grepl("PC", colnames(eig_vec))]
  eig_values[[Sample]] <- data.frame(
    PC = PCs, 
    eigen_value = read_lines(paste0(eig_vectors_dir, eig_val)), 
    sample = unname(sample_names[Sample])) %>%
    mutate(
      eigen_value = as.numeric(eigen_value),
      eigen_percentage = eigen_value / sum(eigen_value)
      )
  
  #get the eigen vector file and store in a list
  eig_vectors[[Sample]] <- eig_vec %>%
    bind_rows(.id = "table") %>%
    select(-`#FID`) %>%
    rename(barcode = IID) %>%
    mutate(barcode = gsub("_", "", barcode)) %>%
    mutate(sample = sample_names[Sample])
}

eig_values <- bind_rows(eig_values)

##merge them in a single data frame and append strain information from cells_meta
eig_vectors <- bind_rows(eig_vectors) %>%
  left_join(
    cells_meta %>%
      select(-barcode) %>%
      rename(barcode = simple_barcode) %>% 
      group_by(barcode) %>%
      mutate(n = n()) %>% #there is one barcode (AGGTCCAACGCTACTACACCATCTTTCCACTC) which is present in both scDNA_AT_01 and scDNA_AT_02 libraries. without the number ID I can't tell which is which, so remove it here
      filter(n == 1) %>%
      select(barcode, sample, strain, cell_or_background),
    by = c("sample", "barcode")
  ) %>%
  mutate(strain = ifelse(is.na(strain), "doublet", strain)) %>% #doublets are not present in cells_meta, so they are assigned as NAs. Fix it here.
  mutate(sample = factor(sample, levels = sample_names[sample_names %in% sample])) #to rearrange sample by the order specified in `sample_names`



#plot the PCA####
##first create axis labels for the PCs
eig_labels <- eig_values %>%
  filter(PC %in% c("PC1", "PC2")) %>%
  group_by(PC) %>%
  summarise(min_percent = min(eigen_percentage), max_percent = max(eigen_percentage)) %>%
  mutate(label = paste0(PC, " (", round(min_percent * 100, 2), "% - ", round(max_percent * 100, 2), "%", ")")) %>%
  select(PC, label) %>%
  deframe()

figure_5A <- eig_vectors %>%
  ggplot(aes(x = PC1, y = PC2, fill = strain))+
  geom_point(size = 2, shape = 21, stroke = 0.1)+
  facet_wrap(vars(sample), scales = 'free', nrow = 1)+
  guides(color = guide_legend(override.aes = list(size = 2)))+
  scale_fill_manual(values = c(strain_colors, doublet = "grey"))+
  labs(x = eig_labels[["PC1"]], y = eig_labels[["PC2"]])+
  theme_bw()


#plot distribution of NA values####
##basically this summarises how many sites could not be resolved in each cell
mat_list <- list()
for(Sample in names(sample_names)){
  #get the allele depth matrix.
  mat <- vroom::vroom(paste0(allele_depth_dir, "/", Sample, "_ad_matrix.tsv"), col_types = cols(.default = col_character())) %>%
    select(-last_col()) %>% # don't know why but the last column is just an empty column, so remove it here. 
    mutate(site = paste0(chromosome, "_", position)) %>%
    select(-chromosome, -position) %>%
    column_to_rownames("site") %>%
    select( #remove doublets
      any_of(
        cells_meta %>%
          filter(strain != "doublet") %>%
          pull(simple_barcode)
        )
      ) 
  
  Sample <- sample_names[Sample]
  mat_list[[Sample]] <- mat
}


#calculate for each site the percentage of cells which could not resolve that site (allele depth = 0,0)
NA_plot <- list()
for(Sample in names(mat_list)){
  NA_plot[[Sample]] <- data.frame(sample = Sample, mean_NA = rowMeans(mat_list[[Sample]] == "0,0", na.rm = TRUE))
}

#convert it to a single data frame
NA_plot <- NA_plot %>%
  bind_rows() %>%
  mutate(sample = factor(sample, levels = sample_names[sample_names %in% sample])) 

#plot it
figure_5B <- NA_plot %>%
  ggplot(aes(y = 1-mean_NA, x = sample))+
  geom_violin(aes(fill = sample))+
  geom_boxplot(width = 0.075, outlier.size = 0.1, size = 0.2)+
  scale_y_continuous(name = "fraction of cells with\nconfident variant call", breaks = c(0:10)/10)+
  scale_x_discrete(name = NULL)+
  scale_fill_manual(values = sample_colors)+
  guides(fill = "none")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#plot heatmap with alleles####
##first we need to get the drug-resistance positions to make sure they are included in the heatmap
#get strain anotation
strains <- cells_meta %>%
  select(simple_barcode, strain) %>%
  deframe()

#get the snp effect data for sample_6. It will be used for both plots
snpeff_data <- read.delim(paste0(drugs_data_dir, "/sample_6.filtered.drugs.snpeff.csv")) %>%
    fuzzy_inner_join( #append the drug markers data
      drug_data,
      by = c("CHROM" = "chrom", "POS" = "start", "POS" = "end"),
      match_fun = list(`==`, `>=`, `<=`)
    ) %>%
  filter(!ANN.0..EFFECT %in% c("synonymous_variant")) %>% #removes synonymous variants
  mutate(row = paste0(CHROM, "_", POS)) #create a column for the row name combining the chromosome and position
  
plot_list <- list()
for(Sample in c("SPC-STD2", "SPC-PTA2")){
  ##get the data frequency matrix
  hm_mat <- mat_list[[Sample]]
  ##make sure drug rows are in the matrix
  rows <- snpeff_data$row #get the rows where drug-resistance allels are
  rows <- c(rows, sample(rownames(hm_mat), size = 1000)) #combine it with a sample of rows (matrix is too large to plot complete)
  rows <- unique(rows) #remove duplicates
  rows <- rows[rows %in% rownames(hm_mat)] #remove rows which are not in the matrix
  hm_mat <- hm_mat[rows,] #subset the matrix.
  
  ##convert it to allele frequency
  ### since this is a complex task I paralellize here
  p <- detectCores() - 1L
  
  res_list <- mclapply(seq_len(ncol(hm_mat)), function(j){
    col <- hm_mat[, j]
    s <- scan(text = col, sep = ",", what = list(integer(), integer()), quiet = TRUE)
    ref <- s[[1]]
    alt <- s[[2]]
    
    freq <- ifelse(ref+alt == 0L, NA_real_, alt / (alt+ref))
  }, 
  mc.cores = p)
  
  res <- do.call(cbind, res_list)
  dimnames(res) <- dimnames(hm_mat)
  hm_mat <- res
  rm(res)
  
  #
  
  #cluster the cells
  hclust_cells <- hclust(daisy(t(hm_mat), metric = "gower"), method = "ward.D2")
  
  #arrange rows by mean allele frequency (placing the the strain defining ones on top)
  vars_order <- rowMeans(hm_mat, na.rm = TRUE)
  vars_order <- sort(vars_order)
  vars_order <- names(vars_order)
  hm_mat <- hm_mat[vars_order,]
  
 
  major_hm <- ggheatmap(hm_mat)+
    scale_x_discrete(name = NULL, breaks = NULL)+
    scale_y_discrete(name = NULL, breaks = NULL)+
    scale_fill_gradient2(name = "allele frequency",
                         breaks = c(0, 0.5, 1), 
                         limits = c(0,1), 
                         low = snp_colors[1], mid = snp_colors[2], high = snp_colors[3], 
                         midpoint = 0.5, 
                         na.value = "white")+
    theme(legend.position = "bottom")+
    anno_top(size = 0.3)+
    align_phylo(as.phylo(hclust_cells))+
    scale_y_continuous(name = NULL, breaks = NULL)+
    patch_titles(top = Sample)+
    anno_top() + 
    ggalign(data = strains[colnames(hm_mat)], size = 0.3)+ 
    geom_tile(aes(y = 1, fill = factor(value))) +
    scale_fill_manual(name = "Strain", values = strain_colors)+ 
    theme_void()+
    theme(legend.position = "top")
  
  if(Sample == "SPC-PTA2"){
    #get the index of the rows that contain the drug-resistance variants, to highlight them with arrows in 
    row_idx <- which(rownames(hm_mat) %in% snpeff_data$row)
    major_hm <- major_hm +
      anno_right(size = 0.1)+
      ggalign(data = seq_len(nrow(hm_mat))) +
      geom_star(data = data.frame(y = row_idx, x = 0), aes(x = x, y = y), inherit.aes = FALSE, starshape = 11, angle = 270, size = 2, fill = "white", colour = "black")+
      scale_x_continuous(name = NULL, breaks = NULL, expand = c(0,0), limits = c(0,1))+
      coord_cartesian(clip = "off")
  }
 
  
  plot_list[[Sample]] <- major_hm
}

figure_5C <- ggalign::align_plots(!!!plot_list, nrow = 1, guides = "tblr") + layout_tags(NULL)

#now make a second heatmap just for the drug-resistance genes only in sample 6
##since sample 6 was the last one to used above, we can reuse the hm_mat from the for loop above.
hm_drug <- hm_mat[row_idx,]
#convert it to factor
hm_drug[hm_drug >= 0.9] <- 2L
hm_drug[hm_drug <= 0.1] <- 0L
hm_drug[hm_drug > 0.1 & hm_drug < 0.9] <- 1L
hm_drug[,] <- apply(hm_drug, 2, factor)

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
  mutate(label = ANN.0..HGVS_P) %>%
  mutate(label = gsub("p.", "", label)) %>%
  dplyr::select(row, label) %>%
  deframe()

rownames(hm_drug) <- gene_labels[rownames(hm_drug)]
names(gene_groups) <- gene_labels[names(gene_groups)]

#chose rows to make the additional bar plot
rows_to_bar <- rowMeans(hm_drug == "1", na.rm = TRUE)
rows_to_bar <- which(rows_to_bar > 0.05)

#prepare a dataframe for the barplot
bar_data <- hm_drug[rows_to_bar, ] %>%
  as_tibble(rownames = "locus") %>%
  pivot_longer(cols = -locus, names_to = "barcode", values_to = "allele") %>%
  mutate(strain = strains[barcode]) %>%
  filter(strain == "HU3") 

drug_hm <- ggheatmap(hm_drug)+
  scale_x_discrete(name = NULL, breaks = NULL)+
  scale_y_continuous(position = "right")+
  scale_fill_manual(name = "Allele", breaks = c(0, 1, 2), labels = c(`0` = "Homozygous reference", `1` = "Heterozygous", `2` = "Homozygous alternative"), values = snp_colors, na.value = "white")+
  theme(axis.text.y = element_text(size = 8),
        strip.text.y.right = element_blank(),
        strip.text.y.left = element_text(angle = 0, hjust = 0, vjust = 0.5, size = 10), #this is needed to keep panel labels (gene groups) as default `theme_ggalign()` removes them. 
        strip.placement = "outside", #allow labels to be draw outside the box
        strip.clip = "off", 
        legend.position = "bottom") + #prevent labels from being clipped. 
  facet_grid(switch = "y") +
  anno_left()+
  align_group(group = gene_groups[rownames(hm_drug)])+
  anno_top(size = 0.3)+
  align_phylo(as.phylo(hclust_cells))+
  scale_y_continuous(name = NULL, breaks = NULL)+
  anno_top() + 
  ggalign(data = strains[colnames(hm_mat)], size = 0.3)+ 
  geom_tile(aes(y = 1, fill = factor(value))) +
  scale_fill_manual(name = "Strain", values = strain_colors)+ 
  theme_void()+
  theme(legend.position = "top")
#anno_right()+
#ggmark(mark_line(rows_to_bar))+
#geom_bar(data = bar_data, aes(fill = allele, x = locus), position = "fill", show.legend = FALSE, na.rm = TRUE)+
#scale_fill_manual(name = "Allele", breaks = c(0, 1, 2), labels = c(`0` = "Homozygous reference", `1` = "Heterozygous", `2` = "Homozygous alternative"), values = snp_colors)+
#theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#ggtitle("Fraction in HU3")

figure_5D <- drug_hm

##now we plot the phylogenetic trees made in sample 6
#get the trees
file_names <- list.files(path = "inputs/trees", pattern = ".*sample_6.*.nwk")
file_names <- file_names[!grepl("both", file_names)] #remove the tree with both strains together
trees <- paste0("inputs/trees/", file_names)
trees <- lapply(trees, read.tree)
names(trees) <- file_names

plot_list <- list()
for(file in names(trees)){
  tree <- trees[[file]]
  p <- ggtree(tree, layout = "ellipse")
  p$data <- cbind(p$data, cells_meta[p$data$label,])
  title <- unique(p$data$strain)
  title <- title[!is.na(title)]
  p <- p +
    geom_tippoint(aes(fill = strain), size = 1.5, shape = 21, stroke = 0.3) +
    geom_rootpoint(color = "black") +
    scale_fill_manual(values = strain_colors)+
    theme_tree2()+
    ggtitle(title)+
    guides(fill = "none", color = "none")+
    coord_cartesian(clip = "off")
  
  if(title == "HU3"){
    arrow_plot <- ggplot() +
      geom_segment(aes(x = 0, xend = 1, y = 0, yend = 0),
                   arrow = arrow(length = unit(0.2, "cm")),
                   linewidth = 0.5) +
      xlim(-0.2, 1.2) +
      ylim(-0.5, 0.5) +
      theme_void()
    
    zoom_p <- viewClade(p, node = 150, xmax_adjust = -1500) + ggtitle(NULL) 
    p <- p + geom_hilight(node = 150, fill = "skyblue", alpha = 0.2)
    p <- ggalign::align_plots(p, arrow_plot, zoom_p, nrow = 1, widths = c(0.6, 0.1, 0.3)) + layout_tags(NULL)
  }
  
  
  
  plot_list[[file]] <- p
}

figure_5E <- ggalign::align_plots(!!!plot_list, nrow = 1, widths = c(0.4, 0.6))
figure_5E

top <- ggalign::align_plots(free_border(figure_5A), figure_5B, nrow = 1, widths = c(0.85, 0.15))
mid <- ggalign::align_plots(figure_5C, figure_5D, nrow = 1, widths = c(0.6, 0.4))
final <- ggalign::align_plots(top, mid, NULL, figure_5E, ncol = 1, heights = c(0.15, 0.7, 0.03, 0.27))
final <- final + layout_tags("A") + layout_theme(plot.tag = element_text(size = 16))


#save the final figure panel####
plot_scale <- 1.8
ggsave("figure_5.pdf", plot = final, width = 8.27 * plot_scale, height = 9 * plot_scale)

#open it
if (Sys.info()["sysname"] == "Darwin") {
system2("open", args = "figure_5.pdf", wait = FALSE)
} else if (Sys.info()["sysname"] == "Linux") {
system2("xdg-open", args = "figure_5.pdf", wait = FALSE)
}