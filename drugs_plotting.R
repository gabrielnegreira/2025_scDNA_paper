library('vcfR')
library('adegenet')
library('pheatmap')
library('stringr')
library('RColorBrewer')
library('reshape2')
library('ggplot2')
library(openxlsx)


## input parameters
src_dir = '/Users/pmonsieurs/programming/leishmania_scDNA_atrandi/results/drugs/'
setwd(src_dir)


## for filtered SNP vcf file
sample = 5
vcf_file = paste0('sample_', sample, '.filtered.drugs.snpeff.vcf')
snpeff_file = paste0('sample_', sample, '.filtered.drugs.snpeff.csv')

drug_gene_file = '/Users/pmonsieurs/programming/leishmania_susl/data/drug_resistance_Ldon_v2.xlsx'



#### 1. overall analysis ####


snp_colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(3)
breaks = c(0,.66,1.33,2.00)

# ---- 1.1 input data ----

## use vcfR to read data
vcf <- read.vcfR(vcf_file, verbose = TRUE)
queryMETA(vcf)
queryMETA(vcf, element = 'ANN')
head(getFIX(vcf))
vcf@gt[1:5,]

## us genlight object
x <- vcfR2genlight(vcf)
gt_matrix = t(as.data.frame(x))
dim(gt_matrix)

## remove those SNPs that are the same (or NA) in all samples
gt_matrix_clean <- gt_matrix[apply(gt_matrix, 1, function(row) {
  unique_values <- unique(na.omit(row)) # Remove NA values and get unique values
  length(unique_values) > 1            # Keep rows with more than one unique value
}), ]

## create heatmap
pheatmap(gt_matrix_clean,
         fontsize_col = 8,
         show_rownames = FALSE,
         cluster_cols = TRUE,
         color = snp_colors,
         breaks = breaks, 
         treeheight_row = FALSE)



# ---- 1.2 focus on SNPs in coding regions ----

## only select those position where the SNP has in impact, which means that it 
## should occur in the coding region, and the SNP should have an impact (i.e. 
## a nonsense or missense mutation). The impact of those SNPs and indels can be
## assessed using the SnpEff approach, 

## read in the info field from the vcf file, which contains the information of 
## SNPeff
info_data = INFO2df(vcf)
ann_data = str_split_fixed(info_data$ANN, "\\|", n=10)
head(ann_data)
table(ann_data[,2])

## select only those SNP with a substantial impact (e.g. ("missense_variant", 
## "stop_gained") or remove the one with a non-harmful deletion (e.g. 
## "downstream_gene_variant", "synonymous_variant", "upstream_gene_variant")
variants_omit = c("start_lost", "downstream_gene_variant", "synonymous_variant", "upstream_gene_variant", "intergenic_region")
vcf_impact = vcf[! ann_data[,2] %in% variants_omit, ]

## create heatmap with impact SNPs
x_impact <- vcfR2genlight(vcf_impact)
gt_impact_matrix = t(as.data.frame(x_impact))


## remove those SNPs that are the same (or NA) in all samples
gt_impact_matrix <- gt_impact_matrix[apply(gt_impact_matrix, 1, function(row) {
  unique_values <- unique(na.omit(row)) # Remove NA values and get unique values
  length(unique_values) > 1            # Keep rows with more than one unique value
}), ]

pheatmap(gt_impact_matrix,
         fontsize_col = 8,
         show_rownames = FALSE,
         show_colnames = FALSE,
         cluster_cols = TRUE,
         color = snp_colors,
         breaks = breaks, 
         treeheight_row = FALSE)



## rename the rownames so that the gene is included
snpeff_data = read.csv(snpeff_file, sep="\t")
drug_data = read.xlsx(drug_gene_file, sheet = "Sheet1")
snp_chroms = sapply(rownames(gt_impact_matrix), function(x) strsplit(x, "_")[[1]][1])
snp_pos = sapply(rownames(gt_impact_matrix), function(x) strsplit(x, "_")[[1]][2])

snp_information = data.frame(chrom = snp_chroms,
                                 pos = snp_pos)


snp_information$drug = NA
snp_information$mutation = NA
snp_information$class = NA
snp_information$effect = NA
snp_information$gene = NA


for (i in 1:nrow(gt_impact_matrix)) {
  snp_info = rownames(gt_impact_matrix)[i]
  row_index = which(rownames(snp_information) == snp_info)
  
  ## find chrom and position in the matrix
  chrom = snp_information$chrom[row_index]
  pos = as.integer(snp_information$pos[row_index])
  
  ## find the drug
  drug_name = drug_data[drug_data$start <= pos & drug_data$end >= pos & chrom == chrom,]$Drug_code
  snp_information[row_index,]$drug = drug_name
  
  gene = drug_data[drug_data$start <= pos & drug_data$end >= pos & chrom == chrom,]$Gene_code
  snp_information[row_index,]$gene = gene
  
  
  
  ## find the mutation
  snp_information[row_index,]$mutation = snpeff_data[snpeff_data$CHROM == chrom & snpeff_data$POS == pos,]$ANN.0..HGVS_P
  snp_information[row_index,]$effect = snpeff_data[snpeff_data$CHROM == chrom & snpeff_data$POS == pos,]$ANN.0..EFFECT
  
}

## create new rownames
rownames(gt_impact_matrix) == rownames(snp_information)
rownames(gt_impact_matrix) =  paste0(snp_information$chrom, "_",
                                     snp_information$pos, "_",
                                     snp_information$class, "_",
                                     snp_information$drug, "_",
                                     snp_information$mutation)

rownames(gt_impact_matrix) =  paste0(snp_information$drug, "_",
                                     snp_information$gene, "_",
                                     snp_information$chrom, "_",
                                     snp_information$pos, "::",
                                     snp_information$effect, "_",
                                     snp_information$mutation)                                         

p = pheatmap(gt_impact_matrix,
         fontsize_col = 6,
         fontsize_row = 10,
         show_rownames = TRUE,
         show_colnames = FALSE,
         cluster_cols = TRUE,
         cluster_rows = FALSE,
         color = snp_colors,
         breaks = breaks, 
         treeheight_row = FALSE,
         treeheight_col = FALSE)
p

png_file = paste0(src_dir, 'snp_heatmap_drug_resistance_sample', sample, '.png')
ggsave(file = png_file, plot=p, width=16, height=9, dpi=300)


## find the frameshift variant in AQP1
head(gt_impact_matrix)
colnames(gt_impact_matrix)[which(gt_impact_matrix[1,] == 2)]



## make the sum of the different types of SNPs per position
for (i in 1:nrow(gt_impact_matrix)) {
  plot_data = as.data.frame(gt_impact_matrix[i,])
  colnames(plot_data) = c('snp_type')
  plot_data$snp_type = gsub("0", '0_homozygote ref', plot_data$snp_type)
  plot_data$snp_type = gsub("1", '1_heterozygote SNP', plot_data$snp_type)
  plot_data$snp_type = gsub("2", '2_homozygote SNP', plot_data$snp_type)
  plot_data$count = 1
  plot_data$sample = rownames(gt_impact_matrix)[i]
  p = ggplot(data = plot_data, aes(x=sample, y=count, fill=snp_type)) + 
    geom_bar(position="stack", stat = "identity") + 
    theme_bw()
  
  png_out = paste0(src_dir, 'sample_', sample, '_snptype_', rownames(gt_impact_matrix)[i], '.png')
  ggsave(plot = p, file = png_out, dpi=300, width=4, height=9)

}

