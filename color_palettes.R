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
names(sample_colors) <- c(paste("Sample", c(1:6)), "10X BPK081", "10X HU3", "10X")

#colors for snp heatmaps
snp_colors <- rev(brewer.pal(7, name = "RdYlBu")[c(1,4,7)])
names(snp_colors) <- c(0, 1, 2)

#use this command to visualize a color palette
colors <- sample_colors
#barplot(rep(1, length(colors)),col = colors,border = NA, space = 0, names.arg = seq_along(colors), las = 1)