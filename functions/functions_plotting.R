# analysis functions

CalcTSS=function(inter, thresh = 0.5)
{
  
  aW=which(inter$value==1 & inter$pred1 >=thresh)
  bW=which(inter$value==0 & inter$pred1 >=thresh) 
  cW=which(inter$value==1 & inter$pred1 <thresh) 
  dW=which(inter$value==0 & inter$pred1 <thresh) 
  
  a=as.numeric(length(aW)) # TRUE POSITIVE
  b=as.numeric(length(bW)) # FALSE POSITIVE (Value = 0, but incorrectly marked as 1)
  c=as.numeric(length(cW)) # FALSE NEGATIVE (Value = 1, but incorrectly marked as 0)
  d=as.numeric(length(dW)) # TRUE NEGATIVE
  
  Py=a/(a+b)*100
  Pn=d/(c+d)*100
  
  TSS=(a*d-b*c)/((a+c)*(b+d))
  
  TSSd=data.frame(a=a,b=b,c=c,d=d,TSS=TSS)
  
  return(TSSd)
}

CalcTSS_factor=function(inter)
{
  
  aW=which(inter$value==1 & inter$pred1 == 1)
  bW=which(inter$value==0 & inter$pred1 == 1) 
  cW=which(inter$value==1 & inter$pred1 == 0) 
  dW=which(inter$value==0 & inter$pred1 == 0) 
  
  a=as.numeric(length(aW)) # TRUE POSITIVE
  b=as.numeric(length(bW)) # FALSE POSITIVE (Value = 0, but incorrectly marked as 1)
  c=as.numeric(length(cW)) # FALSE NEGATIVE (Value = 1, but incorrectly marked as 0)
  d=as.numeric(length(dW)) # TRUE NEGATIVE
  
  Py=a/(a+b)*100
  Pn=d/(c+d)*100
  
  TSS=(a*d-b*c)/((a+c)*(b+d))
  
  TSSd=data.frame(a=a,b=b,c=c,d=d,TSS=TSS)
  
  return(TSSd)
}

calc_fw_metrics <- function(.data = NULL)
{
  
  # simplify graph
  graph <- igraph::simplify(.data, remove.multiple = T, remove.loops = T)
  # prepare properties
  properties <- list("title" = "Iberian Ponds" , "M.units" = "µgC", "N.units" = "indvs/L")
  # prepare trophic links
  trophic.links <- igraph::as_data_frame(graph, what = "both")$edges %>% 
    dplyr::select(resource = to, consumer = from)
  # prepare nodes
  nodes <- igraph::as_data_frame(.data, what = "both")$vertices %>% 
    dplyr::select(node = name, avg.N, avg.M)
  # load the community
  iberian_community <- cheddar::Community(properties=properties, nodes=as.data.frame(nodes), trophic.links = trophic.links)
  # remove computation limits
  options(cheddarMaxQueue=0)
  #chains <- tryCatch(TrophicChains(iberian_community), error=print)
  
  PBFL = cheddar::PreyAveragedTrophicLevel(iberian_community, include.isolated = F)
  
  #FBFL = cheddar::FlowBasedTrophicLevel(iberian_community, weight.by = NULL, include.isolated = FALSE),
  #LTL = cheddar::LongestTrophicLevel(iberian_community, include.isolated = TRUE)
  #ShortestTrophicLevel(iberian_community, include.isolated=FALSE)
  #ShortWeightedTrophicLevel(iberian_community, include.isolated=TRUE)
  #LongWeightedTrophicLevel(iberian_community, include.isolated=TRUE)
  #CATL <- ChainAveragedTrophicLevel(iberian_community, include.isolated=TRUE)
  #TH <- TrophicHeight(iberian_community, include.isolated=FALSE)
  #TLs <- TrophicLevels(iberian_community, weight.by="B.units", include.isolated=TRUE)

  return(PBFL)
  
}

calc_trophic_vulnerability <- function(.data = NULL)
{
  
  # simplify graph
  graph <- igraph::simplify(.data, remove.multiple = T, remove.loops = T)
  # prepare properties
  properties <- list("title" = "Iberian Ponds" , "M.units" = "µgC", "N.units" = "indvs/L")
  # prepare trophic links
  trophic.links <- igraph::as_data_frame(graph, what = "both")$edges %>% 
    dplyr::select(resource = to, consumer = from)
  # prepare nodes
  nodes <- igraph::as_data_frame(.data, what = "both")$vertices %>% 
    dplyr::select(node = name, avg.N, avg.M)
  # load the community
  iberian_community <- cheddar::Community(properties=properties, nodes=as.data.frame(nodes), trophic.links = trophic.links)

  trophic_vulnerability <- cheddar::TrophicVulnerability(iberian_community)
  
  return(trophic_vulnerability)
  
}

calc_trophic_generality <- function(.data = NULL)
{
  
  # simplify graph
  graph <- igraph::simplify(.data, remove.multiple = T, remove.loops = T)
  # prepare properties
  properties <- list("title" = "Iberian Ponds" , "M.units" = "µgC", "N.units" = "indvs/L")
  # prepare trophic links
  trophic.links <- igraph::as_data_frame(graph, what = "both")$edges %>% 
    dplyr::select(resource = to, consumer = from)
  # prepare nodes
  nodes <- igraph::as_data_frame(.data, what = "both")$vertices %>% 
    dplyr::select(node = name, avg.N, avg.M)
  # load the community
  iberian_community <- cheddar::Community(properties=properties, nodes=as.data.frame(nodes), trophic.links = trophic.links)
  
  trophic_generality <- cheddar::TrophicGenerality(iberian_community)
  
  return(trophic_generality)
  
}

## plot functions 

set_colours_by_category <- function(.data = NULL.) {
  
  ## get guilds
  nodes <- igraph::as_data_frame(.data, "both")$vertices
  nodes_groups <- nodes$category
  num_nodes <- length(nodes_groups)
  nodes_trophic_colors <- nodes_groups
  
  # set colours per group
  nodes_trophic_colors[nodes_groups == "phyto_producer"] <- bar_colors[1]
  nodes_trophic_colors[nodes_groups == "heterotrophic microeukaryot_predator"] <- bar_colors[2]
  nodes_trophic_colors[nodes_groups == "heterotrophic microeukaryot_omnivore"] <- bar_colors[3]
  nodes_trophic_colors[nodes_groups == "heterotrophic microeukaryot_herbivore"] <- bar_colors[4]
  nodes_trophic_colors[nodes_groups == "zoo_herbivore"] <- bar_colors[5]
  nodes_trophic_colors[nodes_groups == "zoo_omnivore"] <- bar_colors[6]
  nodes_trophic_colors[nodes_groups == "zoo_predator"] <- bar_colors[7]
  nodes_trophic_colors[nodes_groups == "macro_herbivore"] <- bar_colors[8]
  nodes_trophic_colors[nodes_groups == "macro_detritivore"] <- bar_colors[9]
  nodes_trophic_colors[nodes_groups == "macro_omnivore"] <- bar_colors[10]
  nodes_trophic_colors[nodes_groups == "macro_predator"] <- bar_colors[11]
  
  return(nodes_trophic_colors)
  
}


set_colours_by_old <- function(.data = NULL., by = "trophic") {
  
  if (by == "trophic") {
    
    # get groups
    nodes <- V(.data)
    nodes_groups <- nodes$group
    num_nodes <- length(nodes_groups)
    nodes_trophic_colors <- nodes_groups
    
    # set colours per group
    nodes_trophic_colors[nodes_groups == "macro"] <- "#E6ED18"
    nodes_trophic_colors[nodes_groups == "zoo"] <- "#78D857"
    nodes_trophic_colors[nodes_groups == "phyto"] <- "#006D80"
    
    # assign ranks
    nodes_trophic_rank <- nodes_trophic_colors
    nodes_trophic_rank[nodes_groups == "phyto"] <- 1
    nodes_trophic_rank[nodes_groups == "zoo"] <- 2
    nodes_trophic_rank[nodes_groups == "macro"] <- 3
    nodes_trophic_rank <- as.numeric(nodes_trophic_rank)
    
  }
  
  if (by == "guild") {
    
    # get guilds
    nodes <- V(.data)
    nodes_groups <- nodes$category
    num_nodes <- length(nodes_groups)
    nodes_trophic_colors <- nodes_groups
    
    # set colours per group
    nodes_trophic_colors[nodes_groups == "phyto_producer"] <- bar_colors[1]
    nodes_trophic_colors[nodes_groups == "heterotrophic microeukaryot_predator"] <- bar_colors[2]
    nodes_trophic_colors[nodes_groups == "heterotrophic microeukaryot_omnivore"] <- bar_colors[3]
    nodes_trophic_colors[nodes_groups == "heterotrophic microeukaryot_herbivore"] <- bar_colors[4]
    nodes_trophic_colors[nodes_groups == "zoo_herbivore"] <- bar_colors[5]
    nodes_trophic_colors[nodes_groups == "zoo_omnivore"] <- bar_colors[6]
    nodes_trophic_colors[nodes_groups == "zoo_predator"] <- bar_colors[7]
    nodes_trophic_colors[nodes_groups == "macro_herbivore"] <- bar_colors[8]
    nodes_trophic_colors[nodes_groups == "macro_detritivore"] <- bar_colors[9]
    nodes_trophic_colors[nodes_groups == "macro_omnivore"] <- bar_colors[10]
    nodes_trophic_colors[nodes_groups == "macro_predator"] <- bar_colors[11]
    
    # assign ranks
    nodes_trophic_rank <- nodes_trophic_colors
    nodes_trophic_rank[nodes_groups == "phyto_producer"] <- bar_colors[1]
    nodes_trophic_rank[nodes_groups == "heterotrophic microeukaryot_predator"] <- bar_colors[2]
    nodes_trophic_rank[nodes_groups == "heterotrophic microeukaryot_omnivore"] <- bar_colors[3]
    nodes_trophic_rank[nodes_groups == "heterotrophic microeukaryot_herbivore"] <- bar_colors[4]
    nodes_trophic_rank[nodes_groups == "zoo_herbivore"] <- bar_colors[5]
    nodes_trophic_rank[nodes_groups == "zoo_omnivore"] <- bar_colors[6]
    nodes_trophic_rank[nodes_groups == "zoo_predator"] <- bar_colors[7]
    nodes_trophic_rank[nodes_groups == "macro_herbivore"] <- bar_colors[8]
    nodes_trophic_rank[nodes_groups == "macro_detritivore"] <- bar_colors[9]
    nodes_trophic_rank[nodes_groups == "macro_omnivore"] <- bar_colors[10]
    nodes_trophic_rank[nodes_groups == "macro_predator"] <- bar_colors[11]
    nodes_trophic_rank <- as.numeric(nodes_trophic_rank)
    
  }
  
  return(list(nodes_trophic_colors, nodes_trophic_rank))
  
}

set_edges_by_interaction_type <- function(.data = NULL) {
  
  edge_list <- as_long_data_frame(.data)
  interaction_type_groups <- edge_list$trait_2_from
  num_edges <- length(interaction_type_groups)
  edges_trophic_colors <- interaction_type_groups
  
  # set colours per group
  edges_trophic_colors[interaction_type_groups == "phyto_producer"] <- bar_colors[1]
  edges_trophic_colors[interaction_type_groups == "heterotrophic microeukaryot_predator"] <- bar_colors[2]
  edges_trophic_colors[interaction_type_groups == "heterotrophic microeukaryot_omnivore"] <- bar_colors[3]
  edges_trophic_colors[interaction_type_groups == "heterotrophic microeukaryot_herbivore"] <- bar_colors[4]
  edges_trophic_colors[interaction_type_groups == "zoo_herbivore"] <- bar_colors[5]
  edges_trophic_colors[interaction_type_groups == "zoo_omnivore"] <- bar_colors[6]
  edges_trophic_colors[interaction_type_groups == "zoo_predator"] <- bar_colors[7]
  edges_trophic_colors[interaction_type_groups == "macro_herbivore"] <- bar_colors[8]
  edges_trophic_colors[interaction_type_groups == "macro_detritivore"] <- bar_colors[9]
  edges_trophic_colors[interaction_type_groups == "macro_omnivore"] <- bar_colors[10]
  edges_trophic_colors[interaction_type_groups == "macro_predator"] <- bar_colors[11]
  
  return(list(edges_trophic_colors))
  
}

## plot functions 
set_colours_by <- function(.data = NULL., by = "category") {
  
  if (by == "category") {
    
    # get category
    nodes <- V(.data)
    nodes_groups <- nodes$category
    num_nodes <- length(nodes_groups)
    nodes_trophic_colors <- nodes_groups
    
    # set colours per group
    nodes_trophic_colors[nodes_groups == "phyto_producer"] <- bar_colors[1]
    nodes_trophic_colors[nodes_groups == "heterotrophic microeukaryot_predator"] <- bar_colors[2]
    nodes_trophic_colors[nodes_groups == "heterotrophic microeukaryot_omnivore"] <- bar_colors[3]
    nodes_trophic_colors[nodes_groups == "heterotrophic microeukaryot_herbivore"] <- bar_colors[4]
    nodes_trophic_colors[nodes_groups == "zoo_herbivore"] <- bar_colors[5]
    nodes_trophic_colors[nodes_groups == "zoo_omnivore"] <- bar_colors[6]
    nodes_trophic_colors[nodes_groups == "zoo_predator"] <- bar_colors[7]
    nodes_trophic_colors[nodes_groups == "macro_herbivore"] <- bar_colors[8]
    nodes_trophic_colors[nodes_groups == "macro_detritivore"] <- bar_colors[9]
    nodes_trophic_colors[nodes_groups == "macro_omnivore"] <- bar_colors[10]
    nodes_trophic_colors[nodes_groups == "macro_predator"] <- bar_colors[11]
    
    # assign ranks
    # nodes_trophic_rank <- nodes_trophic_colors
    # nodes_trophic_rank[nodes_groups == "producer"] <- 1
    # nodes_trophic_rank[nodes_groups == "herbivore"] <- 2
    # nodes_trophic_rank[nodes_groups == "omnivore"] <- 3
    # nodes_trophic_rank[nodes_groups == "predator"] <- 4
    # nodes_trophic_rank[nodes_groups == "detritivore"] <- 5
    # nodes_trophic_rank <- as.numeric(nodes_trophic_rank)
  }
  
  
  if (by == "group") {
    
    # get feeding mode
    nodes <- V(.data)
    nodes_groups <- nodes$group
    num_nodes <- length(nodes_groups)
    nodes_trophic_colors <- nodes_groups
    
    # set colours per group
    nodes_trophic_colors[nodes_groups == "phyto"] <- bar_colors[1]
    nodes_trophic_colors[nodes_groups == "zoo"] <- bar_colors[2]
    nodes_trophic_colors[nodes_groups == "macro"] <- bar_colors[3]
    
    # assign ranks
    #nodes_trophic_rank <- nodes_trophic_colors
    #nodes_trophic_rank[nodes_groups == "phyto"] <- 1
    #nodes_trophic_rank[nodes_groups == "zoo"] <- 2
    #nodes_trophic_rank[nodes_groups == "macro"] <- 3
    #nodes_trophic_rank <- as.numeric(nodes_trophic_rank)
  }
  return(list(nodes_trophic_colors))##nodes_trophic_rank))
}

#' prepare_graphs_plotting
#'
#' @param .edges Edge list 
#' @param .traits Traits dataframe
#'
#' @return Edge list
#' @export
#'
#' @examples prune_full_metaweb(edge_list)
prepare_graphs_plotting <- 
  function(.edges = NULL, .traits = NULL, .threshold_value = 0, .column_name = "species_code") {
    
    #.edges <- hybrid_metaweb_interactions
    #.traits <- hybrid_data_trait_phylo
    #.column_name = "unique_id"
    #.threshold_value = 0
    
    .edges <- .edges[.edges$from != .edges$to, ]
    
    #probability threshold
    .edges <- .edges[.edges$weight > .threshold_value, ]
    
    #remove NAs from body mass
    .traits <- .traits[!is.na(.traits$avg.M),]
    
    #change order of the columns
    #.traits<- .traits[,c(2,1,3,4,5,6,7,8,9,10,11,12,13,14,15)]
    
    #set traits to vertexes
    edge_list_mat <- graph_from_data_frame(.edges, directed = T, vertices = NULL) # TODO ZEYNEP: before it was false
    
    # match species in the species list with those in the trait list
    spp_trait_list <- unique(as.vector(.traits$species_code)) # get species codes from both dataframes
    spp_edge_list <- unique(as_ids(V(edge_list_mat)))
    shared_species <- intersect(spp_trait_list, spp_edge_list)
    edge_list_mat <- induced_subgraph(edge_list_mat, shared_species)
    input_traits_subset <- .traits[.traits$species_code %in% shared_species,]
    # step to make sure the order of the traits and nodes match
    order_vertices <- igraph::as_data_frame(edge_list_mat, what = "vertices")[,1]
    input_traits_subset <- input_traits_subset[match(order_vertices, input_traits_subset$species_code),]
    
    if(.column_name == "unique_id") {
      
      #set traits to vertexes
      edge_list_mat <- graph_from_data_frame(.edges, directed = T, vertices = NULL) # TODO ZEYNEP: before it was false
      # match species in the species list with those in the trait list
      spp_trait_list <- unique(as.vector(.traits$unique_id)) # get species codes from both dataframes
      spp_edge_list <- unique(as_ids(V(edge_list_mat)))
      shared_species <- intersect(spp_trait_list, spp_edge_list)
      edge_list_mat <- induced_subgraph(edge_list_mat, shared_species)
      input_traits_subset <- .traits[.traits$unique_id %in% shared_species,]
      # step to make sure the order of the traits and nodes match
      order_vertices <- igraph::as_data_frame(edge_list_mat, what = "vertices")[,1]
      input_traits_subset <- input_traits_subset[match(order_vertices, input_traits_subset$unique_id),]
      
      }
    
    # attach traits information to each vertexes
    edge_list_mat <-
      edge_list_mat %>%
      set_vertex_attr(., c("group"), value = input_traits_subset$group) %>%
      set_vertex_attr(., c("category"), value = input_traits_subset$category) %>%
      set_vertex_attr(., c("avg.M"), value = input_traits_subset$avg.M) %>%
      set_vertex_attr(., c("avg.N"), value = input_traits_subset$avg.N)
    
    return(edge_list_mat)
  }

# this function is deprecated since 25/08
create_igraph_object <- function(all_edges_table, ib_traits) {
  #to remove self links (if from and to are the same species)
  all_edges_table <- all_edges_table[all_edges_table$from != all_edges_table$to, ]
  
  #probability threshold
  threshold_value = 0
  all_edges_table <- all_edges_table[all_edges_table$weight > threshold_value, ]
  
  #remove NAs from body mass
  ib_traits <- ib_traits[!is.na(ib_traits$avg.M),]
  
  #change order of the columns
  #ib_traits<- ib_traits[,c(2,1,3,4,5,6,7,8,9,10,11,12,13,14,15)]
  
  #set traits to vertexes
  edge_list_mat <- graph_from_data_frame(all_edges_table, directed = F, vertices = NULL)
  
  # match species in the species list with those in the trait list
  spp_trait_list <- unique(as.vector(ib_traits$species_code)) # get species codes from both dataframes
  spp_edge_list <- unique(as_ids(V(edge_list_mat)))
  shared_species <- intersect(spp_trait_list, spp_edge_list)
  edge_list_mat <- induced_subgraph(edge_list_mat, shared_species)
  input_traits_subset <- ib_traits[ib_traits$species_code %in% shared_species,]
  
  # step to make sure the order of the traits and nodes match
  order_vertices <- igraph::as_data_frame(edge_list_mat, what = "vertices")[,1]
  input_traits_subset <- input_traits_subset[match(order_vertices, input_traits_subset$species_code),]
  
  # attach traits information to each vertexes
  edge_list_mat <-
    edge_list_mat %>%
    set_vertex_attr(., c("group"), value = input_traits_subset$group)%>%
    set_vertex_attr(., c("category"), value = input_traits_subset$category) %>%
    set_vertex_attr(., c("avg.M"), value = input_traits_subset$avg.M)%>%
    set_vertex_attr(., c("avg.N"), value = input_traits_subset$avg.N)
  
  return(edge_list_mat)
}



plot_interactions <- function(edge_list_mat, title) {
  # set colours and ranks
  node_colours <- set_colours_by(edge_list_mat, by = "category")[[1]]
  #ranks <- set_colours_by(edge_list_mat, by = "category")[[2]]
  
  # set layout
  layout_web <- layout_nicely(edge_list_mat)
  
  layout_web[, 2] <- log(rescale(as.vector(V(edge_list_mat)$avg.M), c(1, 100)))
  #layout_web[, 2] <- log(as.vector(V(edge_list_mat)$body_mass) + 0.1, 1000)
  
  plot.igraph(edge_list_mat,
              asp = 0.8,
              rescale = TRUE, xlim = c(-1, 1), ylim = c(-1, 1),
              mark.groups = T, mark.col = 0, mark.border = 0,
              edge.color = 2, edge.arrow.size = 0, edge.curved = 0.2,
              edge.width = 1,
              vertex.size = 10,
              vertex.color = node_colours,
              vertex.label = "",
              vertex.label.dist = 2,
              layout = layout_web,
              main=title,
              # "\n(N=", graph.properties$N,
              # ", L=", round(graph.properties$Ltot, digits = 2),
              # ", LD=", round(graph.properties$LD, digits = 2),
              # ", C=", round(graph.properties$C, digits = 2), ")",
              # sep = ""
  )
  #dev.off()
}


