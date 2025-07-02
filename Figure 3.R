#clean the environment####
rm(list = ls(envir = parent.frame()))
gc()

#load libraries####
library(tidyverse)
library(ggridges)
library(ggrepel)
library(ggalluvial)
library(patchwork)
library(ggalign)
library(igraph)
library(ggraph)
library(RGraphSpace)
library(GGally)
library(ape)
library(pegas)
library(visNetwork)
source("scDNA_functions_S3.R")
source("color_palettes.R")

#get inputs####
true_cells <- readRDS("inputs/true_cells.rds")
scDNA10X <- readRDS("inputs/scDNA10X.rds")

##set base theme parameters for ggplot####
update_geom_defaults("point", list(size = 0.5))

#bind cells meta####
cells_meta <- lapply(c(true_cells, scDNA10X), function(x)x$metadata$cells_meta)
cells_meta <- lapply(cells_meta, function(x){
  x$sample <- as.character(x$sample) 
  return(x)})
cols = Reduce(intersect, lapply(cells_meta, names))
cells_meta <- lapply(cells_meta, function(x)x[,cols])
cells_meta <- bind_rows(cells_meta)

#plot figures####
##plot_figure 3A####
###creat a color pallete
breaks <- lapply(true_cells, function(x)x$somies$raw_somy_matrix) #get the somies for each sample
breaks <- breaks[c(1,5)] #subset only samples 1 and 5
breaks <- unique(as.integer(unlist(breaks))) #get all possible somy values
breaks <- c(0:max(breaks)) #set the breaks to 0 to the maximum possible somy.
color_palette <- heat_col(breaks) # create the heat color map.

#create the plot list where both plots will be stored
plot_list <- list()
#make the same plot for samples 2 (index 2) and sample 6 (index 5 because sample 4 was removed from true_cells)
for(i in c(2, 5)){

  #get raw somy matrix
  matrix_to_plot <- true_cells[[i]]$somies$raw_somy_matrix
  
  #get cells' strain
  strains <- true_cells[[i]]$metadata$cells_meta %>%
    rownames_to_column("cell") %>%
    select(cell, strain) %>%
    deframe()
  
  #get the sample
  sample <- unique(true_cells[[i]]$metadata$cells_meta$sample)
  
  #plot density
  dens_plot <- ggplot(data.frame(value = as.vector(matrix_to_plot)), aes(x = value))+
    geom_density()+
    scale_x_continuous(name = "Raw Somies", breaks = c(1:100))
  
  #plot core heatmap
  #plot core heatmap
  hm_plot <- ggheatmap(matrix_to_plot) + #start the heatmap
    scale_fill_gradientn(name = "Raw Somies", colors = color_palette, breaks = breaks, limits = c(0, max(breaks))) + #set the fill scale
    scale_x_discrete(labels = NULL, breaks = NULL)+ #remmoves cell labels
    anno_left(size = 0)+ #create an empty left annotation(needed for the next step)
    align_order(rev(rownames(matrix_to_plot)))+ #reverse the order of the plot
    anno_top(size = 0.2) + #create an anotation space in the top
    free_border(dens_plot, borders = "l") + #align the density plot on top of the matrix
    patch_titles(top = paste("Sample", sample))+ #add a title to the top annotation
    anno_top() + #create another annotation space at the top
    align_dendro(method = "ward.D2", size = 0.2)+ #add the dendogram (which also reorders the columns)
    theme_void()+ #remove irrelevant elements of the dendogram plot
    anno_top() + #create a second annotation space, below the dendogram and on top of the heatmap
    ggalign(data = strains[colnames(matrix_to_plot)], size = 0.1)+ #initiate a ggplot object withy the strains as data.
    geom_tile(aes(y = 1, fill = factor(value))) + #add annotation tiles.
    scale_fill_manual(name = "Strain", values = strain_colors)+ #set the annotation colors
    theme_void()
  
  plot_list <- c(plot_list, hm_plot) 
  
  }

#align both figures
figure_3A <- ggalign::align_plots(!!!plot_list, guides = "r")

##plot figures 3B####
#get karyotypes for each strain
net_plot_list <- list()
karyo_plot_list <- list()
for(strain_to_plot in c("BPK081", "HU3")){
  
  #get the cells of that strain and which are not outliers, and not polyploid
  cells <- true_cells[[5]]$metadata$cells_meta %>%
    filter(strain == strain_to_plot) %>% 
    filter(!is.na(karyotype) & !is_outlier) %>%
    filter(scale_factor <= 2.5) %>%
    rownames()
  
  #get the integer somy values for those cells
  karyotypes <- true_cells[[5]]$somies$int_somy_matrix[,cells]
  
  #mape each cell to its respective karyotype
  cells_karyo <- t(karyotypes) %>%
    as.data.frame() %>%
    unite(karyotype, starts_with("Ld"), sep = "_", remove  = FALSE) %>%
    select(karyotype) %>%
    rownames_to_column("cell") %>%
    deframe()
  
  #create a summary table for each karyotype
  karyo_list <- t(karyotypes) %>%
    as.data.frame() %>%
    unite(karyotype, starts_with("Ld"), sep = "_", remove  = FALSE) %>%
    group_by(karyotype) %>%
    summarise(ncells = n())%>%
    mutate(prop = ncells/sum(ncells)) %>%
    arrange(desc(prop)) %>%
    mutate(karyo_id = paste0("kar", row_number()))
  
  #calculate pairwise distances
  karyo_dist_matrix <- karyo_dist(t(karyotypes))
  
  
  #build the network
  network <- mst(karyo_dist_matrix)
  
  #make a nodes list
  nodes <- data.frame(id = attr(network, "label")) %>%
    mutate(karyotype = cells_karyo[id],
           karyo_id = karyo_list$karyo_id[match(karyotype, karyo_list$karyotype)],
           karyo_ncells = karyo_list$ncells[match(karyotype, karyo_list$karyotype)])
  
  #add unique colors to each cell based on their karyotypes
  colors <- nodes %>%
    filter(karyo_ncells > 1) %>%
    pull("karyotype") %>%
    unique() %>%
    create_colors()
  
  nodes$color <- colors[nodes$karyotype]
  nodes <- nodes %>%
    mutate(color = ifelse(is.na(color), "#f0f0f0", color))
  
  #make the edge list
  #we start by connecting cells with identical karyotypes to a single founder
  edges_karyo <- nodes %>%
    mutate(to = id) %>%
    group_by(karyotype) %>%
    mutate(from = sample(id, size = 1)) %>%
    ungroup() %>%
    filter(to != from) %>%  
    transmute(from, to, width = 1, somy_changes = 0) %>%
    mutate(karyo_from = cells_karyo[from]) %>%
    mutate(karyo_to = cells_karyo[to]) 
  
  #now we use the mst to define edges conecting cells from different karyotypes
  edges_cross <- as.data.frame(network[,])
  colnames(edges_cross) <- c("from", "to", "somy_changes")
  edges_cross$from <- attr(network, "label")[edges_cross$from]
  edges_cross$to <- attr(network, "label")[edges_cross$to]
  
  # add karyotype information
  edges_cross <- edges_cross %>%
    mutate(karyo_from = cells_karyo[from]) %>%
    mutate(karyo_to = cells_karyo[to]) %>%
    filter(karyo_from != karyo_to) 
  
  # Combine both edge lists
  edges <- bind_rows(edges_karyo, edges_cross) 
  
  #convert to igraph object
  g <- graph_from_data_frame(edges[,c("from", "to", "somy_changes")], vertices = nodes, directed = FALSE)
  E(g)$weight <- ifelse(E(g)$somy_changes == 0, 5, 10+E(g)$somy_changes)
  E(g)$label <- ifelse(E(g)$somy_changes <= 1, NA, E(g)$somy_changes)
  E(g)$color <- ifelse(E(g)$somy_changes == 0, "lightgrey", 
                       ifelse(E(g)$somy_changes == 1, "black", "orange"))
  
  
  #set cell groups based on karyotypes
  V(g)$group <- ifelse(V(g)$karyo_ncells == 1, NA, V(g)$karyotype)
  
  # set the network layout
  lay <- layout_with_kk(g, coords = layout_with_fr(g), weights = E(g)$weight)
  
  #check the standard igraph plot
  plot(g, layout = lay, vertex.label = NA, vertex.size = 5)
  
  #now plot with ggraph
  
  net_plot_list[[strain_to_plot]] <- ggraph(g, layout = lay)+
    geom_edge_link(aes(label = ifelse(is.na(label), "", label), color = I(color)))+
    geom_node_point(aes(fill = I(color)), size = 1.5, shape = 21)+
    theme_void()+
    ggtitle(strain_to_plot)+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

figure_3B <- cowplot::plot_grid(plotlist = net_plot_list, ncol = 1)

##make figure 3D####
to_plot <- cells_meta %>%
  filter(strain != "doublet") %>%
  filter(!is.na(karyotype)) %>%
  group_by(strain, sample, karyotype) %>%
  summarise(ncells = n()) 

to_plot %>%
  group_by(sample, strain) %>%
  mutate(karyo_position = as.integer(fct_reorder(karyotype, desc(ncells)))) %>%
  mutate(karyo_id = paste0("kar", karyo_position)) %>%
  filter(karyo_position <= 5) %>%
  separate(karyotype, sep = "_", into = rownames(true_cells[[5]]$somies$int_somy_matrix), remove = FALSE) %>%
  pivot_longer(cols = contains("Ld"), names_to = "chromosome", values_to = "somy") %>%
  mutate(sample = factor(sample, levels = c(1,2,3,5,6, "10X"))) %>%
  ggplot(aes(x = karyo_id, y = fct_rev(chromosome), fill = somy))+
  geom_tile(width = 0.98, col = "#2a5686")+
  labs(x = "Top 5 karyotypes in each data set", y = "Chromosome", fill = "Somy")+
  facet_grid(rows = vars(strain), cols = vars(sample))+
  scale_fill_manual(values = heat_col(c(1:15)))

to_plot <- to_plot %>%
  pivot_wider(names_from = c(strain, sample), values_from = ncells) %>%
  pivot_longer(cols = -karyotype, values_to = "ncells") %>%
  mutate(strain = gsub("_.*", "", name)) %>%
  mutate(sample = gsub(".*_", "", name)) %>%
  select(strain, sample, karyotype, ncells)



#make it pairwise
to_plot$`10X` <- NA
for(Strain in unique(to_plot$strain)){
  ncells_in_10X <- to_plot %>%
    filter(strain == Strain & sample == "10X") %>%
    select(karyotype, ncells) %>%
    deframe()
  to_plot <- to_plot %>%
    mutate(`10X` = ifelse(strain == Strain, ncells_in_10X[karyotype], `10X`))
}

  
separate(karyotype, sep = "_", into = rownames(karyotypes), remove = FALSE) %>%
  



plot <- to_plot %>%
  filter(sample != "10X") %>%
  rename(Atrandi = ncells) %>%
  pivot_longer(cols = c(Atrandi, `10X`), names_to = "experiment", values_to = "ncells") %>%
  mutate(ncells = ifelse(is.na(ncells), 0, ncells)) %>%
  mutate(experiment = fct_rev(factor(experiment))) %>%
  group_by(strain,sample, experiment) %>%
  mutate(proportion = ncells/sum(ncells)) %>%
  group_by(sample, strain, karyotype) %>%
  mutate(max_prop = max(proportion, na.rm = TRUE)) %>%
  arrange(strain, sample, desc(max_prop)) %>%
  group_by(strain, sample, experiment) %>%
  mutate(position = row_number()) %>%
  mutate(color = ifelse(position <= 10, karyotype, NA)) %>%
  mutate(sample = paste("Sample", sample)) %>%
  mutate(experiment = c(Atrandi = "Test", `10X` = "10X scCNV")[experiment]) %>%
  ggplot(aes(x = fct_rev(experiment), y = proportion, alluvium = karyotype, stratum = karyotype, group = karyotype, fill = color))+
  geom_alluvium(decreasing = FALSE, color = NA)+
  #geom_flow(decreasing = FALSE, color = NA)+
  #geom_stratum(decreasing = FALSE, alpha = 0.8, color = NA)+
  guides(fill = "none")+
  labs(x = "Experiment")+
  scale_x_discrete(expand = c(0.1, 0.1))+
  scale_y_continuous(expand = c(0,0))+
  facet_grid(rows = vars(strain), cols = vars(sample))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), 
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 18), 
        strip.text = element_text(size = 14))

figure_3D <- plot+scale_fill_manual(values = create_colors(plot$data$color, palette = "colorblind friendly"))

#since figure 3A is not compatible with patchwork, I have to use ggalign to align it to the figure_3B
#however it does not allow adding the A/B tags, so I had to do it manually via Inkscape later.
figure_3_panel <- align_plots(figure_3A , figure_3B, ncol = 1, heights = c(0.8, 0.2), guides = "t")

#if I try via patchwork, alignment is off (but tags are added)
top <- wrap_plots(wrap_elements(ggalignGrob(figure_3A)), figure_3B, nrow = 1, widths = c(0.75, 0.25))
bottom <- figure_3B + figure_3C + plot_layout(nrow = 1, widths = c(0.4, 0.6))
final <- wrap_plots(wrap_elements(ggalignGrob(figure_3A)), bottom, ncol = 1, heights = c(0.7, 0.3))+
  plot_annotation(tag_levels = 'A')&
  theme(plot.tag = element_text(size = 18, face = "bold"),
        plot.tag.position   = c(0, 1),    # top-left corner
        plot.tag.background = element_rect(fill = NA, colour = NA),
        axis.title = element_text(size = 10),
        text = element_text(size = 10))


#save the final figure panel####
plot_scale <- 1.8
ggsave("figure_3.pdf", plot = final, width = 8.27 * plot_scale, height = 9 * plot_scale)

#open it
system2('open', args = "figure_3.pdf", wait = FALSE)

