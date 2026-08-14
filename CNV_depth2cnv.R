library(ggplot2)
library(openxlsx)
library(pheatmap)
library(gridExtra)
library(reshape2)
library(stringr)
library(patchwork)

src_dir = '/Users/pmonsieurs/programming/leishmania_scDNA_atrandi/results/cnv/cnv_per_gene/'
cnv_files = list.files(src_dir, pattern = ".cnv.csv")
out_dir = paste0(src_dir, 'figures/')


## read in the meta data and create an additional column name
## read in meta data. ID_code_short is now a separate column in the 
## metadata excel and does not need to be constructed in the script
meta_data_file = '/Users/pmonsieurs/programming/leishmania_scDNA_atrandi/data/cells_meta.xlsx'
meta_data = read.xlsx(meta_data_file)
head(meta_data)

samples = c(1,2,3,5,6)
cnv_data_list = list()

for (sample in samples) {

  barcodes = meta_data[meta_data$sample == sample, ]$barcode_updated
  first = 1
  
  for (barcode in barcodes) {
    cnv_file = paste0(src_dir, barcode, '.cnv.csv')

    if (! file.exists(cnv_file)) {
      print(paste0("file ", cnv_file, "does not exist"))
      next
    }else{
      print(cnv_file)
    }
    
    cnv_new_data = read.csv(paste0(cnv_file))
    

    
    ## concatenate the new data to the existing data frame
    if (first == 1) {
      cnv_data = cnv_new_data
      first = 0
    }else{
      cnv_data = cbind.data.frame(cnv_data, cnv_new_data$ratio)
    }
    
    ## create the column names
    sample_name = gsub(".cnv.csv", "", cnv_file)
    sample_name = gsub(src_dir, "", sample_name)
    colnames(cnv_data)[length(colnames(cnv_data))] = sample_name
    
  }

  ## clean up the rows with only NA values
  dim(cnv_data)
  head(cnv_data)
  
  cnv_data_list[[sample]] = cnv_data[,5:ncol(cnv_data)]
  
  
  ## filter out those CNV values that are higher than the cutoff
  cutoff_cnv = 10
  
  cnv_data_up = cnv_data[rowSums(cnv_data[,5:ncol(cnv_data)] > cutoff_cnv, na.rm=TRUE) > 0,]
  rownames(cnv_data_up) = cnv_data_up$gene_id
  
  annotation_data = as.data.frame(meta_data[match(colnames(cnv_data_up[,5:ncol(cnv_data_up)]), meta_data$barcode_updated),]$strain)
  rownames(annotation_data) = colnames(cnv_data_up[,5:ncol(cnv_data_up)])
  colnames(annotation_data) = c('strain')
  
  ## create heatmap with all strains combined
  breaks = seq(0,5,0.05)
  
  data_mat <- t(cnv_data_up[, 5:ncol(cnv_data_up)])
  data_mat <- data_mat[complete.cases(data_mat), ]
  
  heat = pheatmap(t(cnv_data_up[,5:ncol(cnv_data_up)]),
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           show_colnames = TRUE,
           fontsize_col = 10,
           fontsize_row = 3,
           treeheight_col = 0,
           treeheight_row = 0,
           annotation_row = annotation_data,
           breaks = breaks)
  
  png_file = paste0(out_dir, 'heatmap_cnv_sample_', sample, '.png')
  ggsave(heat, file = png_file, width=16, height = 9)
  
}


## analysis on the number of zero values for the rebuttal of Gabriel

chrom_medians = data.frame()
nr_genes_zero = data.frame()
for (sample in samples) {
  
  barcodes = meta_data[meta_data$sample == sample, ]$barcode_updated
  first = 1
  
  for (barcode in barcodes) {
    cnv_file = paste0(src_dir, barcode, '.cnv.csv')
    
    if (! file.exists(cnv_file)) {
      print(paste0("file ", cnv_file, "does not exist"))
      next
    }else{
      print(cnv_file)
    }
    
    cnv_new_data = read.csv(paste0(cnv_file))
    
    ## extract median values for the full chromosome
    chrom_medians_new = unique(cnv_new_data[,c('chrom', 'chrom_median')])
    chrom_medians_new$sample = sample
    chrom_medians_new$cell = barcode
    chrom_medians = rbind(chrom_medians, chrom_medians_new)
    
    
    ## extract the number of genes with zero values
    genes_zero = sum(cnv_new_data$median_coverage == 0)
    nr_genes_zero = rbind(nr_genes_zero, c(sample, genes_zero))
    
  }
}   


## do distribution of zero genes
colnames(nr_genes_zero) = c('sample', 'genes_zero')
nr_genes_zero$sample = as.factor(nr_genes_zero$sample)

head(nr_genes_zero)
ggplot(data = nr_genes_zero, aes(x=genes_zero)) + 
  geom_density(aes(group = sample, colour = sample)) + 
  theme_bw()


## do some stats on the median chromosome values and gene values
chrom_medians$sample = as.factor(chrom_medians$sample)
ggplot(data = chrom_medians, aes(x= chrom_median)) + 
  geom_density(aes(group = sample, colour = sample)) + 
  theme_bw()

dim(chrom_medians)/36

