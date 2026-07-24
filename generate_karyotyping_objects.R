#clean the environment
rm(list = ls(envir=parent.frame()))
gc()

#libraries####
library(xlsx)
library(tidyverse) 
source("git_modules/scDNA_tools/scDNA_functions.R")#this contains most of the functions used for the analysis

#get inputs####
counts_lib1 <- read.delim("inputs/atrandi_count_matrices/scDNA_AT_01_merged.txt")
counts_lib2 <- read.delim("inputs/atrandi_count_matrices/scDNA_AT_02_merged.txt")
atrandi_barcodes <- read.xlsx("inputs/atrandi_count_matrices/Atrandi split-pool barcodes.xlsx", sheetIndex = 1)
metadata <- read.xlsx("inputs/atrandi_count_matrices/strain_distinction.xlsx", sheetIndex = 1)
sample_metadata <- read.xlsx("inputs/atrandi_count_matrices/sample metadata.xlsx", sheetIndex = 1) %>%
  separate_rows(`Barcode.A.wells`, sep = ",")
bins_meta <- read.delim("inputs/atrandi_count_matrices/bin_gc_and_mappability_20kb_bins.tsv", header = TRUE)
out_dir <- "inputs/karyotyping_objects/"

#set parameters####
min_reads_threshold <- 20000

#format inputs####
#OBS: to create a scDNA object we need a count matrix, as well as a cell metadata and a bin metadata data frame for each sample.

##remove underscores from barcodes
#we start by setting some standards for cell barcodes. For instance, they should have only ATCG characters, folowed by a numeric character indicating from which library the cell originates.
colnames(counts_lib1)[-c(1:3)] <- gsub("_", "", colnames(counts_lib1)[-c(1:3)])
colnames(counts_lib2)[-c(1:3)] <- gsub("_", "", colnames(counts_lib2)[-c(1:3)])

##append a number to the cell barcode indicating from which library that cell comes from
colnames(counts_lib1)[-c(1:3)] <- paste0(colnames(counts_lib1)[-c(1:3)], "_1") 
colnames(counts_lib2)[-c(1:3)] <- paste0(colnames(counts_lib2)[-c(1:3)], "_2") 

#combine both libraries as a single count matrix####
count_matrix <- cbind(counts_lib1, counts_lib2[,-c(1:3)]) %>%
  filter(chrom != "Ld37") %>%
  group_by(chrom) %>%
  mutate(rowname = paste0(chrom, "_", row_number())) %>%
  column_to_rownames("rowname") %>%
  dplyr::select(-chrom, -start, -end) %>%
  as.matrix()

#create cells metadata####
cells_meta <- metadata %>%
  mutate(cell_barcode = gsub("_", "", cell_barcode)) %>%
  mutate(cell_id = row_number()) %>%
  mutate(cell_barcode = paste0(cell_barcode, "_", as.numeric(factor(library)))) %>%
  mutate(barcode = cell_barcode) %>%
  column_to_rownames("cell_barcode") %>%
  mutate(barcode_A = substr(barcode, 25,32),
         barcode_B = substr(barcode, 17,24),
         barcode_C = substr(barcode, 9,16),
         barcode_D = substr(barcode, 1,8)) %>%
  mutate(well_A = atrandi_barcodes$Position.in.plate[match(barcode_A, atrandi_barcodes$Barcode)],
         well_B = atrandi_barcodes$Position.in.plate[match(barcode_B, atrandi_barcodes$Barcode)],
         well_C = atrandi_barcodes$Position.in.plate[match(barcode_C, atrandi_barcodes$Barcode)],
         well_D = atrandi_barcodes$Position.in.plate[match(barcode_D, atrandi_barcodes$Barcode)],) %>%
  cbind(., sample_metadata[match(.$well_A,sample_metadata$Barcode.A.wells),]) %>%
  rename(sample = Sample) %>%
  dplyr::select(-Barcode.A.wells)

#prepare bins_meta####
bins_meta <- bins_meta %>%
  filter(chrom != "Ld37") %>%
  group_by(chrom) %>%
  mutate(bin = paste0(chrom, "_", row_number()), is_mappable = mappability >= 0.7, rownames = bin) %>%
  column_to_rownames("rownames") %>%
  rename(chromosome = chrom, gc_content = X5_pct_gc) %>%
  ungroup()

#create a list of scDNA objects, one per sample####
scDNAobj_list <- list()
for(sample_to_get in sort(unique(cells_meta$sample))){
  cells <- rownames(filter(cells_meta, sample == sample_to_get))
  scDNAobj_list[[sample_to_get]] <- build_scDNAobj(count_matrix = count_matrix[,cells], cells_meta = cells_meta[cells,], bins_meta = bins_meta)
}

#fix doublet detection threshold
scDNAobj_list <- lapply(scDNAobj_list, function(x){
  x$metadata$cells_meta <- x$metadata$cells_meta %>%
    mutate(strain = ifelse(fraction_HU3 > 0.05 & fraction_HU3 < 0.85, "doublet", strain)) 
  return(x)
})

#add relevant info to metadata
for(i in seq_along(scDNAobj_list)){
  scDNAobj_list[[i]]$metadata$cells_meta <- scDNAobj_list[[i]]$metadata$cells_meta %>%
    mutate(experiment = "Atrandi") 
  
  scDNAobj_list[[i]]$metadata$bins_meta <- scDNAobj_list[[i]]$metadata$bins_meta %>%
    mutate(experiment = "Atrandi",
           sample = names(scDNAobj_list)[i])
}

#perform the analysis in each object####
scDNAobj_list <- scDNAobj_list %>%
  lapply(tag_true_cells, plot = FALSE) %>%
  lapply(correct_counts, vars_to_correct = c("gc_content")) %>%
  lapply(normalize_counts) %>%
  lapply(tag_outlier_bins, additional_vars_to_check = c(gc_content = "upper", mappability = "lower", sim_unique_read_count = "both")) %>%
  lapply(calc_cells_ICF, max_chromo = 5) %>%
  lapply(tag_outlier_cells, vars_to_check = c(n_reads = "lower", ICF_score = "upper", ICCV = "upper"))

#save object to disk
saveRDS(scDNAobj_list, paste0(out_dir, "/all_SPCs.rds"))

#now get the true cells only####
true_cells <- scDNAobj_list[c(1,2,3,5,6)] %>%
  lapply(subset_cells, cell_or_background == "cell") %>%
  lapply(subset_cells, is_outlier == FALSE) %>%
  lapply(subset_cells, strain != "doublet") %>%
  lapply(subsample_reads, target_total = min_reads_threshold) %>%
  lapply(correct_counts, vars_to_correct = c("gc_content")) %>%
  lapply(tag_outlier_bins, additional_vars_to_check = c(gc_content = "upper", mappability = "lower", sim_unique_read_count = "both", mean_corrected_counts = "both")) %>%
  lapply(calc_cells_ICF, max_chromo = 5) %>%
  lapply(tag_outlier_cells, vars_to_check = c(n_reads = "lower", ICF_score = "upper", ICCV = "upper")) %>%
  lapply(calc_somy, int_method = "GMM") %>%
  lapply(summarise_karyotypes) %>%
  lapply(karyo_network)

#save it
saveRDS(true_cells, paste0(out_dir, "/true_cells.rds"))

#now we do the analysis for the 10X datasets####
scDNA10X <- list()
scDNA10X[["10X BPK081"]] <- build_scDNAobj(h5 = "inputs/10Xcellranger_h5_files/BPK081.h5")
scDNA10X[["10X HU3"]] <- build_scDNAobj(h5 = "inputs/10Xcellranger_h5_files/HU3.h5")

#add relevant info to metadata
for(i in seq_along(scDNA10X)){
  scDNA10X[[i]]$metadata$cells_meta <- scDNA10X[[i]]$metadata$cells_meta %>%
    mutate(experiment = "10X",
           sample = names(scDNA10X)[i], 
           strain = gsub("10X ", "", sample)) 
  
  scDNA10X[[i]]$metadata$bins_meta <- scDNA10X[[i]]$metadata$bins_meta %>%
    mutate(experiment = "10X",
           sample = names(scDNA10X)[i], 
           strain = gsub("10X ", "", sample))
}

#remove cells with less than 20.000 reads (these are considered not good for karyotyping)
scDNA10X <- scDNA10X %>%
  lapply(subset_cells, n_reads >= min_reads_threshold)

#now set numeric id in barcodes of HU3 to 2

scDNA10X[["10X HU3"]]$metadata$cells_meta <- scDNA10X[["10X HU3"]]$metadata$cells_meta %>%
  rownames_to_column("rowname") %>%
  mutate(rowname = gsub("-1", "-2", rowname),
         barcode = gsub("-1", "-2", barcode)) %>%
  column_to_rownames("rowname")

colnames(scDNA10X[["10X HU3"]]$counts$raw_counts) <- gsub("-1", "-2", colnames(scDNA10X[["10X HU3"]]$counts$raw_counts))

#now subset the BPK081 dataset so it has the same number of cells as the HU3
cells_to_keep <- scDNA10X[["10X BPK081"]]$metadata$cells_meta$barcode
#sample the cells
cells_to_keep <- sample(cells_to_keep, size = nrow(scDNA10X[["10X HU3"]]$metadata$cells_meta), replace = FALSE)
#subset the object
scDNA10X[["10X BPK081"]] <- subset_cells(scDNA10X[["10X BPK081"]], barcode %in% cells_to_keep)

#now combine both datasets as a single one
counts <- lapply(scDNA10X, function(x) x$counts$raw_counts) %>%
  do.call(cbind, .)
cells_meta <- lapply(scDNA10X, function(x) x$metadata$cells_meta) %>%
  bind_rows()
bins_meta <- scDNA10X$`10X BPK081`$metadata$bins_meta
cells_meta$sample <- "10X"
bins_meta$sample <- "10X"
scDNA10X <- build_scDNAobj(count_matrix = counts, cells_meta = cells_meta, bins_meta = bins_meta)

#run the analysis
scDNA10X <- scDNA10X %>% 
  subsample_reads(target_total = min_reads_threshold) %>% 
  correct_counts(vars_to_correct = c("gc_content")) %>% 
  normalize_counts() %>% 
  tag_outlier_bins(additional_vars_to_check = c(gc_content = "upper", mappability = "lower", mean_corrected_counts = "both")) %>% 
  calc_cells_ICF(max_chromo = 5) %>% 
  tag_outlier_cells(vars_to_check = c(effective_depth_of_coverage = "lower", ICF_score = "upper", ICCV = "upper")) %>% 
  calc_somy(int_method = "GMM") %>% 
  summarise_karyotypes() %>% 
  karyo_network()

#save it to disk
saveRDS(scDNA10X, paste0(out_dir,"/scDNA10X.rds"))
  