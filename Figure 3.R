#clean the environment####
source("clean_environment.R")

#load libraries####
source("git_modules/scDNA_tools/scDNA_functions.R")
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
  #create a dataframe with the karyotypes and the number of cells in each sample
  karyo_df <- cells_meta %>%
    filter(strain == Strain) %>%
    filter(!is.na(karyotype)) %>%
    mutate(sample = factor(sample, levels = sample_names[sample_names %in% sample])) %>% #set the order of the samples to the one specified by the global variable `sample_names`
    group_by(strain, experiment, sample, karyotype) %>%
    summarise(ncells = n()) %>%
    group_by(sample, strain) %>%
    arrange(desc(ncells), .by_group = TRUE) %>%
    mutate(proportion = ncells/sum(ncells)) %>%
    mutate(karyo_position = as.integer(fct_reorder(karyotype, desc(ncells)))) %>%
    mutate(karyo_id = paste0("kar", karyo_position)) %>%
    mutate(karyo_name = paste0(sample, "_", karyo_position)) %>% #create a unique name for each karyotype by joining the sample name with the karyotype position.
    ungroup()
  
  
  #now append the corresponding position of the karyotype in the 10X data
  karyo_df <- karyo_df %>%
    left_join(
      karyo_df %>% #takes the same data frame (karyo_df)
        filter(experiment == "10X") %>% #keepts only the 10X data
        rename(position_in_10X = karyo_position, ncells_in_10X = ncells) %>% #renames the karyotype position to `position_in_10X`
        select(karyotype, position_in_10X, ncells_in_10X), #keep only the relevant information: the karyotype and its position in the 10X data
      by = "karyotype" #join the two dataframes by the `karyotype` column
      ) %>%
    mutate(plot_label = ifelse(is.na(position_in_10X), "not in 10X", paste0("10X_kar.", position_in_10X))) %>% #makes a label indicating the position of the karyotype in the 10X data (used later for plotting) %>%
    mutate(plot_label_color =  ifelse(plot_label %in% paste0("10X_kar.", c(1:5)), "black", 
                                 ifelse(plot_label == "not in 10X", "darkred", "lightgrey")))

  
  #create a matrix for the heatmap (hm) with the copy numbers
  heat_matrix <- karyo_df %>%
    filter(karyo_position <= 5) %>%
    separate(karyotype, sep = "_", into = rownames(true_cells[[5]]$somies$int_somy_matrix), remove = FALSE) %>%
    column_to_rownames("karyo_name") %>%
    select(contains("Ld")) %>%
    t() %>%
    as.matrix()
  
  #convert the values in the matrix to integers
  row_names <- rownames(heat_matrix) #apply removes rownames, so first we save them 
  heat_matrix <- apply(heat_matrix, 2, as.integer)
  rownames(heat_matrix) <- row_names
  rm(row_names)
  
  #reverse the order of the rows (so the plot looks nice)
  heat_matrix <- heat_matrix[rev(rownames(heat_matrix)),]
  
  #get a vector containing all possible somy values in the matrix (used to map the values to colors in the `geom_fill_manual` bellow)
  breaks <- sort(unique(as.vector(heat_matrix)))
  
  #create a named vector for the x-axis labels and their colors
  karyo_labels <- karyo_df %>%
    select(karyo_name, plot_label) %>%
    deframe()
  
  label_colors <- karyo_df %>%
    select(karyo_name, plot_label_color) %>%
    deframe()
  
  #match the labels and colors to the colnames of the matrix
  karyo_labels <- karyo_labels[colnames(heat_matrix)]
  label_colors <- label_colors[colnames(heat_matrix)]
  
  #create a named vector specifying how to group the columns in the matrix (used later in `align_group` function)
  col_groups <- karyo_df %>%
    filter(karyo_name %in% colnames(heat_matrix)) %>%
    select(karyo_name, sample) %>%
    deframe()
  
  col_groups <- col_groups[colnames(heat_matrix)]
  
  #plot the heat map
  heat_map <- ggheatmap(heat_matrix, filling = NULL) +
    geom_tile(aes(fill = factor(value)), width = 0.98, color = "#2a5686") +
    scale_x_continuous(expand = c(0,0), labels = as.character(karyo_labels))+ 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color = label_colors))+
    labs(y = "Chromosome", fill = "Somy")+
    scale_fill_manual(values = heat_col(breaks), breaks = breaks)+
    facet_grid(switch = "x", drop = TRUE)+
    theme(axis.text = element_text(), strip.text.x = element_text(), strip.placement = "outside") + #this fixes the postion of the panel labels to be below column labels, not above it 
    theme(axis.text = element_text(), strip.text = element_text())+
    anno_bottom()+
    align_group(col_groups)+
    anno_top() &
    theme(panel.spacing = unit(10, "pt"))
  
  
  
  #set a color scale for the bar plot to highlight the top 5 karyotypes in the 10X dataset
  karyo_10X_colors <- setNames(
    rev(brewer.pal(5, "Spectral")),
    paste0("10X_kar.", 1:5)
  )
  
  bar_plot <- karyo_df %>%
    filter(karyo_name %in% colnames(heat_matrix)) %>%
    mutate(color_label = ifelse(position_in_10X <= 5, paste0("10X_kar.", position_in_10X), NA)) %>%
    ggplot(aes(x = karyo_position, y = proportion, label = ncells, fill = color_label))+
    #geom_col(color = "black", alpha = 0.5)+
    geom_col(color = "black", fill = "black")+
    geom_text(vjust = -0.3, size = 3.7, color = "white")+
    geom_text(vjust = -0.3, size = 3.5, color = "black")+
    facet_wrap(vars(sample), nrow = 1)+
    coord_cartesian(clip = "off")+
    ggtitle(Strain)+
    scale_x_discrete(name = NULL, breaks = NULL,  expand = c(0,0))+
    #scale_fill_manual(values = karyo_10X_colors, breaks = names(karyo_10X_colors))+
    theme_minimal()+
    theme(panel.grid = element_blank(),
          plot.margin = margin(0, 0, 0, 0),
          legend.key.height = unit(0.1, "cm"),
          strip.text = element_blank(),
          strip.clip = "off", 
          panel.spacing = unit(10, "pt"))
  
  heat_map <- ggalign::align_plots(bar_plot, heat_map, ncol = 1, heights = c(0.2, 0.8), guides = "r") + layout_tags(NULL)
  
  #now we format the data for the flow plot
  #the flow plot requires that for every group, all barcodes are listed, even if they were not present there. In that case they should have a 0 count value.
  #The easiest way for me to do that was to pivot the data wider (which adds the NAs to missing counts) and pivot it back to long, converting NAs to 0.
  flow_plot_data <- karyo_df %>%
    select(karyotype, strain, sample, ncells) %>% #select only the needed columns
    mutate(name = paste(strain, sample, sep = "<br>")) %>% #stores the strain, experiment, and sample information in a single column `name` to help `pivot_wider` later
    select(karyotype, name, ncells) %>% #select just the needed columns
    pivot_wider(names_from = name, values_from = ncells) %>% #this ensures that all karyotypes are mentioned in all samples, with NA values when they were not present.
    pivot_longer(cols = -karyotype, values_to = "ncells") %>% #convert back to long format, but now with NA added to karyotypes when they were missing from a given sample.
    separate(name, into = c("strain", "sample"), sep = "<br>") %>% #reconstruct the `experiment`, `strain`, and `sample` columns stored in `name`
    filter(sample != "10X") %>% #since we don't want to compare 10X datasets against itself
    left_join(                             #this left_join block re-appends the number of cells for a given karyotype in the corresponding the 10X data
      karyo_df %>%
        filter(experiment == "10X") %>%
        select(karyotype, ncells_in_10X),
      by = "karyotype"
    ) %>%
    rename(`SPC-scDNA` = ncells, `10X-scDNA` = ncells_in_10X) %>% #we now want to have a single column named `experiment` specifying if it is from 10X of Atrandi datasets. So first we rename these two columns...
    pivot_longer(cols = c(`SPC-scDNA`, `10X-scDNA`), names_to = "experiment", values_to = "ncells") %>% #then we pivot them wider. Now column `Atrandi` has the count in atrandi data, column `10X` has the count in the 10X data
    mutate(ncells = ifelse(is.na(ncells), 0, ncells)) %>% #convert NA counts to 0 (because the karyotype was found in 0 cells in that sample).
    group_by(strain, sample, experiment) %>%
    mutate(proportion = ncells/sum(ncells)) %>% #convert counts to relative proportions
    ungroup() %>%
    mutate(sample = factor(sample, levels = sample_names[sample_names %in% unique(sample)])) %>% #rearrange the samples according to the order set by the `sample_names` environment variable.
    left_join(
      karyo_df %>%
        filter(sample == "10X" & karyo_position <= 5) %>%
        arrange(desc(ncells)) %>%
        mutate(color_label = paste0("10X_kar.", karyo_position)) %>%
        select(karyotype, color_label),
      by = "karyotype"
    )
  
  #finally, make the flow plot  
  flow_plot <- flow_plot_data %>%
    ggplot(aes(x = experiment, y = proportion, alluvium = karyotype, stratum = karyotype, group = karyotype, fill = color_label))+
    geom_alluvium(decreasing = FALSE, color = NA)+
    #guides(fill = "none")+
    scale_x_discrete(name = NULL, expand = c(0.1, 0.1))+
    scale_y_continuous(expand = c(0,0))+
    facet_wrap(vars(sample), ncol = 1)+
    theme(legend.position = "right")+
    scale_fill_manual(values = karyo_10X_colors, breaks = names(karyo_10X_colors))+
    labs(fill = "Top 5 karyotypes\nin 10X dataset")
  
  figures[[length(figures) + 1]] <- ggalign::align_plots(heat_map, free_border(flow_plot), nrow = 1, widths = c(0.8, 0.2))
}

final <- ggalign::align_plots(figures[[1]], figures[[2]], figures[[3]], ncol = 1, heights = c(0.3, 0.35, 0.35))
final <- final + layout_tags("A") + layout_theme(plot.tag = element_text(size = 16))

#save the final figure panel####
plot_scale <- 1.8
ggsave("figure_3.pdf", plot = final, width = 8.27 * plot_scale, height = 9.5 * plot_scale)

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
