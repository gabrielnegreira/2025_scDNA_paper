#clean the environment####
source("clean_environment.R")

#load libraries####
source("scDNA_functions_S3.R")
source("color_palettes.R")
library(tidyverse)
library(ggridges)
library(ggrepel)
library(ggalluvial)
library(ggalign)
library(igraph)
library(ggraph)
library(RGraphSpace)
library(GGally)
library(ape)
library(pegas)
library(visNetwork)

#get inputs####
true_cells <- readRDS("inputs/karyotyping_objects/true_cells.rds")
scDNA10X <- readRDS("inputs/karyotyping_objects/scDNA10X.rds")
cells_meta <- read_delim("inputs/cell_qc/cells_meta.tsv") %>%
  column_to_rownames("rowname") %>%
  filter(cell_or_background == "cell" & strain != "doublet" & sample != "Sample 4") %>%
  mutate(sample = factor(sample_names[sample], levels = sample_names))

names(true_cells) <- sample_names[names(true_cells)]

##set base theme parameters for ggplot####
update_geom_defaults("point", list(size = 0.5))

#plot figures####
figures <- list()
##plot_figure 3A####
###creat a color pallete
breaks <- lapply(true_cells, function(x)x$somies$raw_somy_matrix) #get the somies for each sample
breaks <- breaks[c("SPC-STD2","SPC-PTA2")] #subset only samples 2 and 6
breaks <- unique(as.integer(unlist(breaks))) #get all possible somy values
breaks <- c(0:max(breaks)) #set the breaks to 0 to the maximum possible somy.
color_palette <- heat_col(breaks) # create the heat color map.

#create the plot list where both plots will be stored
plot_list <- list()
#make the same plot for samples 2 (index 2) and SPC-PTA2 (index 5 because sample 4 was removed from true_cells)
for(Sample in c("SPC-STD2", "SPC-PTA2")){

  #get raw somy matrix
  matrix_to_plot <- true_cells[[Sample]]$somies$raw_somy_matrix
  
  #make it 
  #get cells' strain
  strains <- true_cells[[Sample]]$metadata$cells_meta %>%
    rownames_to_column("cell") %>%
    select(cell, strain) %>%
    deframe()
  
  #plot density
  dens_plot <- ggplot(data.frame(value = as.vector(matrix_to_plot)), aes(x = value))+
    geom_density()+
    scale_x_continuous(name = "Raw Somies", breaks = c(1:100))+
    scale_y_continuous(name = NULL, breaks = NULL)
  
  #plot core heatmap
  hm_plot <- ggheatmap(matrix_to_plot, filling = NULL) + #start the heatmap
    geom_tile(aes(fill = value))+
    scale_fill_gradientn(name = "Raw Somies", colors = color_palette, breaks = breaks, limits = c(0, max(breaks))) + #set the fill scale
    scale_x_discrete(labels = NULL, breaks = NULL)+ #remmoves cell labels
    anno_left(size = 0)+ #create an empty left annotation(needed for the next step)
    align_order(rev(rownames(matrix_to_plot)))+ #reverse the order of the plot
    anno_top(size = 0.3) + #create an anotation space in the top
    free_border(dens_plot, borders = "l") + #align the density plot on top of the matrix
    patch_titles(top = Sample)+ #add a title to the top annotation
    anno_top() + #create another annotation space at the top
    align_dendro(method = "ward.D2", size = 0.3)+ #add the dendogram (which also reorders the columns)
    theme_void()+ #remove irrelevant elements of the dendogram plot
    anno_top() + #create a second annotation space, below the dendogram and on top of the heatmap
    ggalign(data = strains[colnames(matrix_to_plot)], size = 0.3)+ #initiate a ggplot object withy the strains as data.
    geom_tile(aes(y = 1, fill = factor(value))) + #add annotation tiles.
    scale_fill_manual(name = "Strain", values = strain_colors)+ #set the annotation colors
    theme_void()
  
  plot_list <- c(plot_list, hm_plot) 
}

#align both figures
figures[[length(figures) + 1]] <- ggalign::align_plots(!!!plot_list, guides = "r")


plot_list <-list()
for(Strain in c("BPK081", "HU3")){
  
  ##plot figure 3B####
  to_plot <- cells_meta %>%
    filter(strain == Strain) %>%
    filter(!is.na(karyotype)) %>%
    group_by(strain, experiment, sample, karyotype) %>%
    summarise(ncells = n()) %>%
    group_by(sample, strain) %>%
    arrange(desc(ncells), .by_group = TRUE) %>%
    mutate(proportion = ncells/sum(ncells)) %>%
    mutate(karyo_position = as.integer(fct_reorder(karyotype, desc(ncells)))) %>%
    mutate(karyo_id = paste0("kar", karyo_position)) %>%
    ungroup()
  
  karyo_in_10X <- to_plot %>%
    filter(experiment == "10X") %>%
    select(karyotype, karyo_position) %>%
    deframe()
  
  df <- to_plot %>%
    filter(karyo_position <= 5) %>%
    mutate(position_in_10X = karyo_in_10X[karyotype]) %>%
    mutate(position_in_10X = factor(ifelse(is.na(position_in_10X), "not present",
                                           ifelse(position_in_10X > 5, "not in top 5", position_in_10X)))) %>%
    mutate(name = paste0(sample, "_", karyo_position))
  
  hm <- df %>%
    separate(karyotype, sep = "_", into = rownames(true_cells[[5]]$somies$int_somy_matrix), remove = FALSE) %>%
    column_to_rownames("name") %>%
    select(contains("Ld")) %>%
    t() %>%
    as.matrix()
  
  
  #split the matrix by sample
  col_groups <- df %>%
    select(name, sample) %>%
    deframe()
  
  sample_order <- as.character(sample_names)
  
  # groups present in your plot (from the columns of hm)
  present <- unique(col_groups[colnames(hm)])
  sample_order_present <- sample_order[sample_order %in% present]
  
  # ordered factor **only with present levels**
  grp <- factor(col_groups[colnames(hm)], levels = sample_order_present)
 
  
  #reverse the order of the rows
  hm <- hm[rev(rownames(hm)),]
  breaks <- sort(unique(as.vector(hm)))
  heat_map <- ggheatmap(hm, filling = NULL) +
  geom_tile(aes(fill = value), width = 0.98, color = "#2a5686") +
  scale_x_continuous(expand = c(0,0), labels = rep(paste0("kar.", c(1:5)), times = 6))+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(y = "Chromosome", fill = "Somy")+
  scale_fill_manual(values = heat_col(breaks), breaks = breaks)+
  facet_grid(switch = "x", drop = TRUE)+
  theme(axis.text = element_text(), strip.text.x = element_text(), strip.placement = "outside") + #this fixes the postion of the panel labels to be below column labels, not above it 
  theme(axis.text = element_text(), strip.text = element_text())+
  anno_bottom()+
  align_group(grp)+
  anno_top() &
  theme(panel.spacing = unit(10, "pt"))
  
  bar_plot <- df %>%
    ggplot(aes(x = karyo_position, y = proportion, label = ncells, fill = position_in_10X))+
    geom_col(color = "black")+
    geom_text(vjust = -0.3, size = 3.7, color = "white")+
    geom_text(vjust = -0.3, size = 3.5, color = "black")+
    facet_wrap(vars(sample), nrow = 1)+
    #scale_y_continuous(limits = c(0, 1))+
    coord_cartesian(clip = "off")+
    #{if(Strain == "BPK081")guides(fill = "none")}+
    ggtitle(Strain)+
    scale_x_discrete(name = NULL, breaks = NULL,  expand = c(0,0))+
    scale_fill_manual(values = c(setNames(scico(9, palette = "lajolla")[3:7], c(1:5)), `not present` = "white", `not in top 5` = "lightgrey"))+
    theme_minimal()+
    theme(panel.grid = element_blank(),
          plot.margin = margin(0, 0, 0, 0),
          legend.key.height = unit(0.1, "cm"),
          strip.text = element_blank(),
          strip.clip = "off", 
          panel.spacing = unit(10, "pt"))
  
  hm <- ggalign::align_plots(bar_plot, heat_map, ncol = 1, heights = c(0.2, 0.8), guides = "r") + layout_tags(NULL)
  
  
  ##plot the flow plot####
  to_plot <- to_plot %>%
    mutate(name = paste(experiment, strain, sample, sep = "<br>")) %>%
    select(-proportion, -karyo_position, -karyo_id, -experiment, -strain, -sample) %>%
    pivot_wider(names_from = name, values_from = ncells) %>%
    pivot_longer(cols = -karyotype, values_to = "ncells") %>%
    separate(name, into = c("experiment", "strain", "sample"), sep = "<br>")
  
  #make it pairwise
  to_plot$`10X` <- NA
  ncells_in_10X <- to_plot %>%
    filter(strain == Strain & experiment == "10X") %>%
    select(karyotype, ncells) %>%
    deframe()
  
  to_plot <- to_plot %>%
    mutate(`10X` = ifelse(strain == Strain, ncells_in_10X[karyotype], `10X`))

  flow_plot <- to_plot %>%
    mutate(sample = factor(sample, levels = sample_names)) %>%
    filter(experiment != "10X") %>%
    select(-experiment) %>%
    rename(Atrandi = ncells) %>%
    pivot_longer(cols = c(Atrandi, `10X`), names_to = "experiment", values_to = "ncells") %>%
    mutate(ncells = ifelse(is.na(ncells), 0, ncells)) %>%
    #mutate(experiment = fct_rev(factor(experiment))) %>%
    group_by(sample, experiment) %>%
    mutate(proportion = ncells/sum(ncells)) %>%
    group_by(sample, strain, karyotype) %>%
    mutate(max_prop = max(proportion, na.rm = TRUE)) %>%
    arrange(strain, sample, desc(max_prop)) %>%
    group_by(strain, sample, experiment) %>%
    mutate(position = row_number()) %>%
    mutate(color = ifelse(position <= 10, karyotype, NA)) %>%
    mutate(experiment = c(Atrandi = "SPC-scDNA", `10X` = "10X-scDNA")[experiment]) %>%
    select(experiment, proportion, karyotype, color, strain) %>%
    ggplot(aes(x = experiment, y = proportion, alluvium = karyotype, stratum = karyotype, group = karyotype, fill = color))+
    geom_alluvium(decreasing = FALSE, color = NA)+
    #geom_flow(decreasing = FALSE, color = NA)+
    #geom_stratum(decreasing = FALSE, alpha = 0.8, color = NA)+
    guides(fill = "none")+
    scale_x_discrete(name = NULL, expand = c(0.1, 0.1))+
    scale_y_continuous(expand = c(0,0))+
    facet_wrap(vars(sample), ncol = 1)+
    theme(legend.position = "right")
  
  flow_plot <- flow_plot + scale_fill_manual(values = create_colors(flow_plot$data$color, palette = "colorblind friendly"))
  #flow_plot <- ggalign::align_plots(NULL, flow_plot, ncol = 1, heights = c(0.2, 0.8)) + layout_tags(NULL)
  
  figures[[length(figures) + 1]] <- ggalign::align_plots(hm, free_border(flow_plot), nrow = 1, widths = c(0.8, 0.2))
}

final <- ggalign::align_plots(figures[[1]], figures[[2]], figures[[3]], ncol = 1, heights = c(0.3, 0.35, 0.35))
final <- final + layout_tags("A") + layout_theme(plot.tag = element_text(size = 16))

#save the final figure panel####
plot_scale <- 1.8
ggsave("figure_3.pdf", plot = final, width = 8.27 * plot_scale, height = 10.5 * plot_scale)

#open it
if (Sys.info()["sysname"] == "Darwin") {
system2("open", args = "figure_3.pdf", wait = FALSE)
} else if (Sys.info()["sysname"] == "Linux") {
system2("xdg-open", args = "figure_3.pdf", wait = FALSE)
}


##plot figure 3D####
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

figures[[length(figures) + 1]] <- ggalign::align_plots(!!!net_plot_list, nrow = 1) + layout_tags(NULL)
