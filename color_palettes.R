library(RColorBrewer)
library(scico) #installed with devtools::install_github("thomasp85/scico")

#set name/value pair to replace sample names for final plots
#this also is used to set factor level to order the samples on the plot, so order here matters!
sample_names <- c(`10X BPK081` = "10X BPK081",
                  `10X HU3` = "10X HU3",
                  `10X` = "10X",
                  `Sample 1` = "SPC-STD1",
                  `Sample 2` = "SPC-STD2",
                  `Sample 3` = "SPC-STD3",
                  `Sample 4` = "SPC-STD4",
                  `Sample 5` = "SPC-PTA1",
                  `Sample 6` = "SPC-PTA2")

#create colors
#this function takes a vector as input and returns a named vector with colors as values and name as elements from the input vector
create_colors <-function(x, 
                         palette = c("default", "all colors", "colorblind friendly", "fancy light", 
                                     "fancy dark", "shades", "tarnish", "pastel", "pimp", 
                                     "intense", "fluo", "red roses", "ochre sand", "yellow lime",
                                     "green mint", "ice cube", "blue ocean", "indigo night", "purple wine"), 
                         seed.use = 123){
  #get required libraries  
  require("hues")
  #test inputs
  palette <- match.arg(palette)
  #create pallete parameters
  palette_param <- list(default = c(hmin = 0, hmax = 360, cmin = 30, cmax = 80, lmin = 35, lmax = 80),
                        `all colors` = c(hmin = 0, hmax = 360, cmin = 0, cmax = 100, lmin = 0, lmax = 100),
                        `colorblind friendly` = c(hmin = 0, hmax = 360, cmin = 40, cmax = 70, lmin = 15, lmax = 85),
                        `fancy light` = c(hmin = 0, hmax = 360, cmin = 15, cmax = 40, lmin = 70, lmax = 100),
                        `fancy dark` = c(hmin = 0, hmax = 360, cmin = 8, cmax = 40, lmin = 7, lmax = 40),
                        shades = c(hmin = 0, hmax = 240, cmin = 0, cmax = 15, lmin = 0, lmax = 100),
                        tarnish = c(hmin = 0, hmax = 360, cmin = 0, cmax = 15, lmin = 30, lmax = 70),
                        pastel = c(hmin = 0, hmax = 360, cmin = 0, cmax = 30, lmin = 70, lmax = 100),
                        pimp = c(hmin = 0, hmax = 360, cmin = 30, cmax = 100, lmin = 25, lmax = 70),
                        intense = c(hmin = 0, hmax = 360, cmin = 20, cmax = 100, lmin = 15, lmax = 80),
                        fluo = c(hmin = 0, hmax = 300, cmin = 35, cmax = 100, lmin = 75, lmax = 100),
                        `red roses` = c(hmin = 330, hmax = 20, cmin = 10, cmax = 100, lmin = 35, lmax = 100),
                        `ochre sand` = c(hmin = 20, hmax = 60, cmin = 20, cmax = 50, lmin = 35, lmax = 100),
                        `yellow lime` = c(hmin = 60, hmax = 90, cmin = 10, cmax = 100, lmin = 35, lmax = 100),
                        `green mint` = c(hmin = 90, hmax = 150, cmin = 10, cmax = 100, lmin = 35, lmax = 100),
                        `ice cube` = c(hmin = 150, hmax = 200, cmin = 0, cmax = 100, lmin = 35, lmax = 100),
                        `blue ocean` = c(hmin = 220, hmax = 260, cmin = 8, cmax = 80, lmin = 0, lmax = 50),
                        `indigo night` = c(hmin = 260, hmax = 290, cmin = 40, cmax = 100, lmin = 35, lmax = 100),
                        `purple wine` = c(hmin = 290, hmax = 330, cmin = 0, cmax = 100, lmin = 0, lmax = 40))
  
  x <- sort(unique(x))
  set.seed(seed.use)
  colors <- hues::iwanthue(n = length(x), 
                           hmin = palette_param[[palette]]["hmin"], 
                           hmax = palette_param[[palette]]["hmax"], 
                           lmin = palette_param[[palette]]["lmin"], 
                           lmax = palette_param[[palette]]["lmax"], 
                           cmin = palette_param[[palette]]["cmin"], 
                           cmax = palette_param[[palette]]["cmax"])
  x <- setNames(colors, x)
  return(x)
}



#colors for strains
strain_colors <- scico(n = 5, palette = "bam")[c(2,4)]
strain_colors <- c(strain_colors, "grey")
names(strain_colors) <- c("BPK081", "HU3", "doublet")

#colors for samples
sample_colors <- c("#f9ce84", "#c48519", "#754a00", "#2d1d00", "#b7d6a3", "#307800", "#6baee5", "#0a4677", "#0a4677")
names(sample_colors) <- c("SPC-STD1", "SPC-STD2", "SPC-STD3", "SPC-STD4", "SPC-PTA1", "SPC-PTA2", "10X HU3", "10X BPK081", "10X")

#colors for snp heatmaps
snp_colors <- rev(brewer.pal(7, name = "RdYlBu")[c(1,4,7)])
names(snp_colors) <- c(0, 1, 2)

#use this command to visualize a color palette
#colors <- sample_colors
#barplot(rep(1, length(colors)),col = colors,border = NA, space = 0, names.arg = seq_along(colors), las = 1)