library(RColorBrewer)
library(scico) #installed with devtools::install_github("thomasp85/scico")

#create a function to map somy values to colors
#heat_col####
#this is a function used to generate the colors for the heatmaps
heat_col <- function(breaks){
  require(colorspace)
  #set the color palette
  colors <- c("#001221", "#002342", "#002342", 
              "#014175", "#035ba3", "#00c3ff", 
              "#00ffee", "#33ff00", "#ccff00", 
              "#fffa00","#ffa600", "#D73027", 
              "#A50026", "#541b1b", "#4d0600")
  
  #get the maximum value
  breaks <- as.integer(breaks)
  n <- max(breaks)
  #by generating 9 colors with this palete, we make sure that values between 0 and 8 are mapped to discrete colors.
  colors <- colorRampPalette(colors)(9)
  #if it asks for more than 9 colors, it will generate the additional colors by darkening the last color
  if(n >= 9){
    for(i in c(10:(n+1))){
      colors[i] <- darken(colors[i-1], amount = 0.4)
    }
  } 
  #now covert values to indices. value 0 should be mapped to the first color, and so on.
  breaks <- breaks+1
  #map the colors to the values
  colors <- colors[breaks]
  names(colors) <- breaks-1
  return(colors)
}


#colors for strains
strain_colors <- scico(n = 5, palette = "bam")[c(2,4)]
strain_colors <- c(strain_colors, "grey")
names(strain_colors) <- c("BPK081", "HU3", "doublet")

#colors for samples
sample_colors <- c("#f9ce84", "#c48519", "#754a00", "#2d1d00", "#b7d6a3", "#307800", "#6baee5", "#0a4677", "#0a4677")
names(sample_colors) <- c(paste("Sample", c(1:6)), "10X BPK081", "10X HU3", "10X")

#colors for snp heatmaps
snp_colors <- rev(brewer.pal(7, name = "RdYlBu")[c(1,4,7)])
names(snp_colors) <- c(0, 1, 2)

#use this command to visualize a color palette
colors <- sample_colors
#barplot(rep(1, length(colors)),col = colors,border = NA, space = 0, names.arg = seq_along(colors), las = 1)