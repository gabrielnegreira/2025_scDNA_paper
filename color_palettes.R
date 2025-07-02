source("scDNA_functions_S3.R")
library(RColorBrewer)
library(scico) #installed with devtools::install_github("thomasp85/scico")


strain_colors <- scico(n = 5, palette = "bam")[c(2,4)]
strain_colors <- c(strain_colors, "grey")
names(strain_colors) <- c("BPK081", "HU3", "doublet")

sample_colors <- scico(n = 8)
names(sample_colors) <- c(1:6, "10X BPK081", "10X HU3")

snp_colors <- brewer.pal(7, name = "RdYlBu")[c(1,4,7)]

names(snp_colors) <- c(0, 1, 2)

#use this command to visualize a color palette
#barplot(rep(1, length(colors)),col = colors,border = NA, space = 0, names.arg = seq_along(colors), las = 1)