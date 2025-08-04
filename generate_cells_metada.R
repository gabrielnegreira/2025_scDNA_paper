#small script to make the full cell metadata
#clean the environment
rm(list = ls(envir=parent.frame()))
gc()

#libraries
library(tidyverse)
library(UpSetR)
library(xlsx)
library(ineq)
source("scDNA_functions_S3.R")

#get the karyotyping object
all_SPCs <- readRDS("inputs/karyotyping_objects/all_SPCs.rds")
true_cells <- readRDS("inputs/karyotyping_objects/true_cells.rds")
scDNA10X <- readRDS("inputs/karyotyping_objects/scDNA10X.rds")
all_SPCs[["10X"]] <- scDNA10X
true_cells[["10X"]] <- scDNA10X
rm(scDNA10X)

#get the additional qc files
qc_files <- list.files("inputs/cell_qc/", pattern = "coverages_combined")
names <- qc_files
qc_files <- paste0("inputs/cell_qc/", qc_files)
qc_files <- lapply(qc_files, read.csv, header = FALSE)
names(qc_files) <- names
rm(names)

#fix the column names in the qc files
qc_files[grep("coverages_combined.csv", names(qc_files))] <- lapply(qc_files[grep("coverages_combined.csv", names(qc_files))], function(x){
  colnames(x) <- c('barcode', 'fraction_1', 'fraction_5', 'fraction_10', 'fraction_25', 'median_coverage', 'mean_coverage', 'min_coverage', 'max_coverage')
  return(x)
})
colnames(qc_files[["coverages_combined_subsampling.csv"]]) <- c('barcode', 'fraction_1_sub', 'fraction_5_sub', 'fraction_10_sub', 'fraction_25_sub', 'median_coverage_sub', 'mean_coverage_sub', 'min_coverage_sub', 'max_coverage_sub')


#remove underscores from barcodes and move them to row names
##our barcodes have a number indicating fromw which library they come. So first I need to add it to the barcodes in the additional files
barcodes <- lapply(all_SPCs, function(x)x$metadata$cells_meta) %>%
  bind_rows() %>%
  rownames_to_column("correct_barcode") %>%
  mutate(simple_barcode = gsub("_.*", "", correct_barcode)) %>%
  mutate(simple_barcode = gsub("-.*", "", correct_barcode)) %>%
  select(simple_barcode, correct_barcode) %>%
  deframe()

#now correct the barcodes in the qc files
qc_files <- lapply(qc_files, function(x){
  x <- x %>%
    rownames_to_column("rownames") %>%
    mutate(barcode = gsub("_|subsampled", "", barcode)) %>%
    filter(barcode %in% names(barcodes)) %>%
    mutate(barcode_correct = barcodes[barcode]) %>%
    mutate(rownames = barcode_correct) %>%
    column_to_rownames("rownames")
  return(x)
})


#calculate gini and mapd
true_cells <- lapply(true_cells, function(x){
  counts <- x$counts$corrected_counts #get corrected counts
  counts <- counts[which(rowSums(counts) > 0),] #remove non mapped regions (especially relevant for 10X as there are artificial Ns in the reference genome there)
  
  #normalize by chromosome
  bins_meta <- x$metadata$bins_meta

  for(chromo in unique(bins_meta$chromosome)){
    bins <- bins_meta %>%
      filter(chromosome == chromo) %>%
      rownames()
    bins <- bins[which(bins %in% rownames(counts))]

    counts[bins, ] <- t(t(counts[bins, ,drop = FALSE]) / colMeans(counts[bins, ,drop = FALSE]))
  }
  
  #normalize by cell
  counts <- apply(counts, 2, normalize, method = "mean")
  
  #calculate gini
  gini <- apply(counts, 2, Gini)
  
  #calculate MAPD
  counts <- counts + 0.1
  counts <- log2(counts)
  mapd <- apply(counts, 2, function(x){
    x <- abs(x[2:length(x)] - x[1:(length(x)-1)])
    x <- median(x)
    return(x)
  })
  x$metadata$cells_meta$gini <- gini[rownames(x$metadata$cells_meta)]
  x$metadata$cells_meta$mapd <- mapd[rownames(x$metadata$cells_meta)]
  return(x)
})

#bind the metadata
cells_meta <- lapply(true_cells, function(x)x$metadata$cells_meta)
cells_meta <- bind_rows(cells_meta)

#add the metadata for the non-true cells
#non_true_meta <- lapply(all_SPCs, function(x)x$metadata$cells_meta)
#non_true_meta <- bind_rows(non_true_meta)
#cells_meta[rownames(non_true_meta ),colnames(non_true_meta )] <- non_true_meta 

#fix some inconsistencies
cells_meta <- cells_meta %>%
  mutate(original_n_reads = ifelse(is.na(original_n_reads), n_reads, original_n_reads))
  
#update some values for 10X samples
cells_meta <- cells_meta %>%
  mutate(cell_or_background = ifelse(experiment == "10X", "cell", cell_or_background))

#bind the dataframes by the barcode column
qc_df <- reduce(qc_files, full_join, by = "barcode_correct") %>%
  column_to_rownames("barcode_correct") 

#bind it to the cells_meta
cells_meta <- cbind(cells_meta, qc_df[rownames(cells_meta),])

#export the cells_meta
cells_meta %>%
  rownames_to_column("rowname") %>%
  write_delim(file = "inputs/cell_qc/cells_meta.tsv")
#write.xlsx(cells_meta, file = "inputs/cell_qc/cells_meta.xlsx")

#visualize gini and mapd
cells_meta %>%
  filter(sample != "Sample 4") %>%
  filter(cell_or_background == "cell") %>%
  filter(!is_outlier) %>%
  ggplot(aes(x = mapd, y = gini, color = strain))+
  geom_point()+
  scale_y_continuous(limits = c(0,1), breaks = c(0:10)/10)+
  facet_grid(cols = vars(sample), scales = "free")