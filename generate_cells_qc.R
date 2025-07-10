#small script to combine the additional cell qc parameters
#libraries
library(tidyverse)
library(UpSetR)
library(xlsx)
library(ineq)

#get the karyotyping object
all_cells <- readRDS("inputs/karyotyping_objects/all_cells.rds")
true_cells <- readRDS("inputs/karyotyping_objects/true_cells.rds")
scDNA10X <- readRDS("inputs/karyotyping_objects/scDNA10X.rds")
true_cells <- c(true_cells, scDNA10X)

#get the additional qc files
qc_files <- list.files("inputs/cell_qc/", pattern = "coverages_combined")
names <- qc_files
qc_files <- paste0("inputs/cell_qc/", qc_files)
qc_files <- lapply(qc_files, read.csv, header = FALSE)
names(qc_files) <- names
rm(names)

#fix the column names in the qc files
colnames(qc_files[["coverages_combined.csv"]]) <- c('barcode', 'fraction_1', 'fraction_5', 'fraction_10', 'fraction_25', 'median_coverage', 'mean_coverage', 'min_coverage', 'max_coverage')
colnames(qc_files[["coverages_combined_subsampling.csv"]]) <- c('barcode', 'fraction_1_sub', 'fraction_5_sub', 'fraction_10_sub', 'fraction_25_sub', 'median_coverage_sub', 'mean_coverage_sub', 'min_coverage_sub', 'max_coverage_sub')


#remove underscores from barcodes and move them to row names
##our barcodes have a number indicating fromw which library they come. So first I need to add it to the barcodes in the additional files
barcodes <- lapply(c(all_cells, scDNA10X), function(x)x$metadata$cells_meta) %>%
  bind_rows() %>%
  rownames_to_column("correct_barcode") %>%
  mutate(simple_barcode = gsub("_.*", "", correct_barcode)) %>%
  select(simple_barcode, correct_barcode) %>%
  deframe()
rm(scDNAobj)

qc_files <- lapply(qc_files, function(x){
  x <- x %>%
    rownames_to_column("rownames") %>%
    mutate(barcode = gsub("_|subsampled", "", barcode)) %>%
    mutate(barcode_correct = barcodes[barcode]) %>%
    mutate(rownames = barcode_correct) %>%
    column_to_rownames("rownames")
  return(x)
})

#recalculate gini and mapd to make sure it is correct
true_cells <- lapply(true_cells, function(x){
  counts <- x$counts$normalized_counts
  gini <- apply(counts, 2, Gini)
  mapd <- apply(counts, 2, function(x){
    x <- abs(x[2:length(x)] - x[1:(length(x)-1)])
    x <- median(x)
    return(x)
  })
  x$metadata$cells_meta$gini_new <- gini[rownames(x$metadata$cells_meta)]
  x$metadata$cells_meta$mapd_new <- mapd[rownames(x$metadata$cells_meta)]
  return(x)
})

#replace old gini and mapd values by new ones (temporary, should remove old gini and mapd values from the scrip)
true_cells <- lapply(true_cells, function(x){
  x$metadata$cells_meta <- x$metadata$cells_meta %>%
    rename(gini = gini_new, mapd = mapd_new)
  return(x)
})

#bind the metadata
cells_meta <- lapply(true_cells, function(x)x$metadata$cells_meta)
cells_meta <- bind_rows(cells_meta)

#visualize if row names are the same across datasets
r_names <- lapply(qc_files, rownames)
r_names <- c(r_names, list(true_cells = rownames(cells_meta)))
upset(fromList(r_names))

#bind the dataframes by the barcode column
qc_df <- reduce(qc_files, full_join, by = "barcode_correct") %>%
  column_to_rownames("barcode_correct") 

#bind it to the cells_meta
cells_meta <- cbind(cells_meta, qc_df[rownames(cells_meta),])


#export the cells_meta
write.xlsx(cells_meta, file = "inputs/cell_qc/cells_meta.xlsx")

