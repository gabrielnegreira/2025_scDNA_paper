#small script to combine the additional cell qc parameters
#libraries
library(tidyverse)
library(UpSetR)
library(xlsx)

#get the karyotyping object
true_cells <- readRDS("inputs/karyotyping_objects/true_cells.rds")
scDNA10X <- readRDS("inputs/karyotyping_objects/scDNA10X.rds")
true_cells <- c(true_cells, scDNA10X)
rm(scDNA10X)

#get the additional qc files
qc_files <- list.files("inputs/cell_qc/", pattern = ".csv")
names <- qc_files
qc_files <- paste0("inputs/cell_qc/", qc_files)
qc_files <- lapply(qc_files, read.csv, header = FALSE)
names(qc_files) <- names
rm(names)


#fix the column names in the qc files
colnames(qc_files[["coverages_combined.csv"]]) <- c('barcode', 'fraction_1', 'fraction_5', 'fraction_10', 'fraction_25', 'median_coverage', 'mean_coverage', 'min_coverage', 'max_coverage')
colnames(qc_files[["coverages_combined_subsampling.csv"]]) <- c('barcode', 'fraction_1_sub', 'fraction_5_sub', 'fraction_10_sub', 'fraction_25_sub', 'median_coverage_sub', 'mean_coverage_sub', 'min_coverage_sub', 'max_coverage_sub')
colnames(qc_files[["mapd.csv"]]) <- c("barcode", "mapd", "sample")
colnames(qc_files[["gini.csv"]]) <- c("barcode", "gini", "sample")
qc_files[["mapd.csv"]] <- qc_files[["mapd.csv"]][-1,]
qc_files[["gini.csv"]] <- qc_files[["gini.csv"]][-1,]

#remove sample columns from mapd and gini dataframes
qc_files[grep("mapd|gini", names(qc_files))] <- lapply(qc_files[grep("mapd|gini", names(qc_files))], function(x){
  x$sample <- NULL
  return(x)
  })

#remove underscores from barcodes and move them to row names
qc_files <- lapply(qc_files, function(x){
  x <- x %>%
    rownames_to_column("rownames") %>%
    mutate(barcode = gsub("_|subsampled", "", barcode)) %>%
    mutate(rownames = barcode) %>%
    column_to_rownames("rownames")
  return(x)
})

#bind the metadata
cells_meta <- lapply(true_cells, function(x)x$metadata$cells_meta)
cells_meta <- bind_rows(cells_meta)
rownames(cells_meta) <- gsub("_.*", "", rownames(cells_meta))

#visualize if row names are the same across datasets
r_names <- lapply(qc_files, rownames)
r_names <- c(r_names, list(true_cells = rownames(cells_meta)))
upset(fromList(r_names))

#bind the dataframes by the barcode column
qc_df <- reduce(qc_files, full_join, by = "barcode") %>%
  column_to_rownames("barcode") 

#bind it to the cells_meta
cells_meta <- cbind(cells_meta, qc_df[rownames(cells_meta),])

#fix numeric values being set to character
cells_meta$gini <- as.numeric(cells_meta$gini)
cells_meta$mapd <- as.numeric(cells_meta$mapd)

#update cells_meta with right name for samples
cells_meta <- cells_meta %>%
  mutate(sample = ifelse(sample == "10X", sample, paste("Sample", sample)))

#export the cells_meta
write.xlsx(cells_meta, file = "inputs/cell_qc/cells_meta.xlsx")

