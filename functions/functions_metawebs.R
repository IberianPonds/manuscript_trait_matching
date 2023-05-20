# metaweb_functions
join_sp <- function(sp, traits) {
  # create a table with group traits
  trait_tab <- data.frame(matrix(0,
                                 nrow = nrow(sp),
                                 ncol = length(levels(traits))
  ))
  for (i in 1:length(levels(traits))) {
    trait_tab[, i] <- rowSums(data.frame(sp[, traits == levels(traits)[i]]))
  }
  colnames(trait_tab) <- levels(traits)
  rownames(trait_tab) <- rownames(sp)
  return(trait_tab)
}

#' Convert vector to presence absences.
#'
#' @param .data Species by site dataframe
#' @param input_traits Traits dataframe
#'
#' @return Edge list
#' @export
#'
#' @examples prune_full_metaweb(edge_list)
pres_abs <- function(x)
{
  if (x > 0) x <- 1 else x <- 0
}

#' prepare_igraph_objects_pond
#'
#' @param .data Species by site dataframe
#' @param input_traits Traits dataframe
#'
#' @return Edge list
#' @export
#'
#' @examples prune_full_metaweb(edge_list)
prepare_igraph_objects_pond <- function(.data = NULL, input_traits = NULL)
{
  
  if (!is.data.frame(.data)) {
    stop("Missing community dataframe")
  }
  
  if (is.null(input_traits)) {
    stop("Missing traits dataframe ")
  }
  
  # test object
  #.data <- pond_matrix_list[[1]]
  input_sample_names <- rownames(.data)
  input_species_names <- colnames(.data[colSums(.data)!=0])
  #input_traits <- species_traits_table
  
  #### get all possible links with pond
  pond_full_metaweb <-
    expand.grid(A = input_species_names, B = input_species_names) %>%
    graph_from_data_frame() %>%
    igraph::simplify(remove.loops = TRUE, remove.multiple = TRUE) # needs to define namespace
  
  # match nodes with traits dataframe and subset graph object
  pond_full_metaweb_nodes <- as.character(as_ids(V(pond_full_metaweb)))
  pond_full_metaweb_traits <- input_traits[input_traits$unique_id %in% pond_full_metaweb_nodes,]
  
  pond_full_metaweb <- induced_subgraph(pond_full_metaweb, !is.na(match(pond_full_metaweb_nodes, pond_full_metaweb_nodes)))
  
  # set node attributes
  pond_full_metaweb <-
    pond_full_metaweb %>%
    set_vertex_attr(name = "group", value = pond_full_metaweb_traits$group) %>%
    set_vertex_attr(name = "avg.N", value = pond_full_metaweb_traits$avg.N) %>%
    set_vertex_attr(name = "avg.M", value = pond_full_metaweb_traits$avg.M) %>%
    set_vertex_attr(name = "category", value = pond_full_metaweb_traits$category)
  
  set_vertex_attr(pond_full_metaweb, name = "group", value = pond_full_metaweb_traits$group)
  
  objectOutputList <- list("graph" = pond_full_metaweb)
  return(objectOutputList)
}

fix_missing_nodes <- function(.data)
{
  
  # get bodymasses
  avgMs <- igraph::as_data_frame(.data, what = "vertices")$avg.M
  
  # identify incomplete nodes
  incomplete_nodes <- V(.data)[is.na(avgMs)]
  
  # remove incomplete nodes
  output_fixed_graph <- .data-incomplete_nodes
  
  return(output_fixed_graph)
  
}

subset_pond_graph <- function(.comm_data = NULL, .metaweb = NULL) {
  
  # species present in pond
  .species_codes <- colnames(.comm_data[as.vector(.comm_data>0)])
  
  # match pond nodes with traits dataframe and subset graph object
  shared_nodes <- intersect(as_ids(V(.metaweb)), .species_codes)
  
  # subset graph object using species present
  subset_metaweb_graph_pond <- induced_subgraph(graph = .metaweb,
                                                vids =  shared_nodes
  )
  
  return(subset_metaweb_graph_pond)
  
}


# function to add detritivores
add_detritus_to_graph <- function(.graph) {
  
  .vertices <- igraph::as_data_frame(.graph, what = "vertices")
  .edges <- igraph::as_data_frame(.graph, what = "edges")
  
  # select detritivores
  .detritivores <- .vertices %>% 
    dplyr::filter(category == "zoo_omnivore" | category == "macro_detritivore") %>% 
    dplyr::select(name) %>% unlist() %>% as.vector()
  
  # create edges for all detritivores
  .edges_plusdetritus <- rbind(.edges, expand_grid(from = .detritivores, to = "detritus", weight = 1))
  
  # create edges for all detritivores
  .vertices_plusdetritus <- rbind(.vertices, data.frame(name = "detritus", group = "detritus", category = "detritus", avg.M = min(.vertices$avg.M), avg.N = 1))
  
  corrected_graph <- graph_from_data_frame(d = .edges_plusdetritus, vertices = .vertices_plusdetritus)
  
  return(corrected_graph)
  
}




#' prepare_igraph_objects
#'
#' @param .data Species by site dataframe
#' @param input_traits Traits dataframe
#' @param .data .data, input_traits, threshold_value = 0
#'
#' @return igraph objects
#' @export
#'
#' @examples prune_full_metaweb(edge_list)
prepare_igraph_objects <- function(.data = NULL, input_traits = NULL, threshold_value = 0)
{
  
  if (!is.data.frame(.data)) {
    stop("Missing community dataframe")
  }
  
  if (is.null(input_traits)) {
    stop("Missing traits dataframe ")
  }
  
  # test objects
  #threshold_value <- 0
  #.data <- species_matrix_MR
  #input_traits <- species_traits_table_per_regions[[1]]
  input_sample_names <- rownames(.data)
  input_species_names <- colnames(.data[colSums(.data)!=0])
  
  # convert to presence absence
  .data[.data > 0] <- 1
  
  # generate co-occurence matrix with crossproduct
  co_mat <- t(.data) %*% as.matrix(.data)
  
  # set diagonal values to 0
  diag(co_mat) <- 0
  
  # create vector with species names species
  spp_names <- rownames(co_mat)
  
  # probably of co-occurence based on the total of potential sites
  max_sites <- dim(.data)[1]
  
  # probability of co-occurence of two species
  prob_co_mat <- co_mat / max_sites
  
  # 4 create igraph object OUTPUT
  edge_list_mat <- graph.adjacency(prob_co_mat, weighted = TRUE)
  
  ## get all links OUTPUT
  all_edges_table <- get.data.frame(edge_list_mat, what = "edges")
  
  # probability threshold
  all_edges_table <- all_edges_table[all_edges_table$weight > threshold_value, ]
  edge_list_mat <- graph_from_data_frame(all_edges_table, directed = T, vertices = NULL)# before directed was False
  
  # match species in the species list with those in the trait list
  spp_trait_list <- unique(as.vector(input_traits$unique_id)) # get species codes from both dataframes
  spp_edge_list <- unique(as_ids(V(edge_list_mat)))
  shared_species <- intersect(spp_trait_list, spp_edge_list)
  edge_list_mat <- induced_subgraph(edge_list_mat, shared_species)
  input_traits_subset <- input_traits[input_traits$unique_id %in% shared_species,]
  
  # step to make sure the order of the traits and nodes match
  order_vertices <- igraph::as_data_frame(edge_list_mat, what = "vertices")[,1]
  input_traits_subset <- input_traits_subset[match(order_vertices, input_traits_subset$unique_id),]
  
  # attach traits information to each vertexes (body sizes)
  edge_list_mat <-
    edge_list_mat %>%
    set_vertex_attr(., c("avg.N"), value = input_traits_subset$avg.N) %>%
    set_vertex_attr(., c("avg.M"), value = input_traits_subset$avg.M) %>%
    set_vertex_attr(., c("group"), value = input_traits_subset$group) %>%
    set_vertex_attr(., c("category"), value = input_traits_subset$category)
  
  # Remove unconnected nodes
  # V(edge_list_mat)$comp <- components(edge_list_mat)$membership
  # edge_list_mat <- induced_subgraph(edge_list_mat, V(edge_list_mat)$comp == 1)
  # graph <- simplify(graph, remove.loops = T, remove.multiple = T)
  
  objectOutputList <- list(
    #"cooccurences" = co_mat,
    "probabilities" = prob_co_mat,
    "species_names" = spp_names,
    "nodes" = V(edge_list_mat),
    "links" = all_edges_table,
    "graph" = edge_list_mat
  )
  return(objectOutputList)
}

#' prepare_igraph_objects_ponds
#'
#' @param .data Species by site dataframe
#' @param input_traits Traits dataframe
#'
#' @return Pruned metaweb
#' @export
#'
#' @examples prune_full_metaweb(edge_list)
prepare_igraph_objects_ponds <- function(.data = NULL, input_traits = NULL, threshold_value = 0)
{
  
  if (!is.data.frame(.data)) {
    stop("Missing community dataframe")
  }
  
  if (is.null(input_traits)) {
    stop("Missing traits dataframe ")
  }
  
  # test objects
  #threshold_value <- 0
  #.data <- pond_matrix_list[[1]]
  #input_traits <- species_traits_table
  input_sample_names <- rownames(.data)
  input_species_names <- colnames(.data[colSums(.data)!=0])
  
  # convert to presence absence
  .data[.data > 0] <- 1

  # generate co-occurence matrix with crossproduct
  co_mat <- t(.data) %*% as.matrix(.data)
  
  # set diagonal values to 0
  diag(co_mat) <- 0
  
  # create vector with species names species
  spp_names <- rownames(co_mat)

  # probably of co-occurence based on the total of potential sites
  max_sites <- dim(.data)[1]
  
  # probability of co-occurence of two species
  prob_co_mat <- co_mat / max_sites
  
  # 4 create igraph object OUTPUT
  edge_list_mat <- graph.adjacency(prob_co_mat, weighted = TRUE)
  
  ## get all links OUTPUT
  all_edges_table <- get.data.frame(edge_list_mat, what = "edges")
  
  # probability threshold
  all_edges_table <- all_edges_table[all_edges_table$weight > threshold_value, ]
  
  edge_list_mat <- graph_from_data_frame(all_edges_table, directed = F, vertices = NULL)
  
  # match species in the species list with those in the trait list
  spp_trait_list <- unique(as.vector(input_traits$unique_id)) # get species codes from both dataframes
  spp_edge_list <- unique(as_ids(V(edge_list_mat)))
  
  shared_species <- intersect(spp_trait_list, spp_edge_list)
  edge_list_mat <- induced_subgraph(edge_list_mat, shared_species)
  
  input_traits_subset <- input_traits[input_traits$unique_id %in% shared_species,]
  
  # step to make sure the order of the traits and nodes match
  order_vertices <- igraph::as_data_frame(edge_list_mat, what = "vertices")[,1]
  input_traits_subset <- input_traits_subset[match(order_vertices, input_traits_subset$unique_id),]
  
  # attach traits information to each vertexes (body sizes)
  edge_list_mat <- set_vertex_attr(edge_list_mat, c("avg.N"), value = input_traits_subset$avg.N)
  edge_list_mat <- set_vertex_attr(edge_list_mat, c("avg.M"), value = input_traits_subset$avg.M)
  edge_list_mat <- set_vertex_attr(edge_list_mat, c("group"), value = input_traits_subset$group)
  edge_list_mat <- set_vertex_attr(edge_list_mat, c("category"), value = input_traits_subset$category)
  
  # Remove unconnected nodes
  # V(edge_list_mat)$comp <- components(edge_list_mat)$membership
  # edge_list_mat <- induced_subgraph(edge_list_mat, V(edge_list_mat)$comp == 1)
  # graph <- simplify(graph, remove.loops = T, remove.multiple = T)
  
  objectOutputList <- list(
    #"cooccurences" = co_mat,
    "probabilities" = prob_co_mat,
    "species_names" = spp_names,
    "nodes" = V(edge_list_mat),
    "links" = all_edges_table,
    "graph" = edge_list_mat
  )
  
  return(objectOutputList)
}

#' Prepare an edge_list for pruning
#' 
#' This function uses the edge list and known traits to 
#' create a new edgelist for each trait combination. This new edgelist
#' is the basis for the pruning functions.
#' 
#' @param .data igraph object to be transformed into an edge list.
#'
#' @return Pruned metaweb
#' @export
#'
#' @examples prune_full_metaweb(edge_list)
prepare_edge_list_for_pruning <- function(.data = NULL)
{
  #.data<- full_metaweb
  if (!is_igraph(.data)) {
    stop("Missing igraph object")
  }
  
  # extract edge list, node list and species list
  selected_edge_list <- igraph::as_data_frame(.data, what = "both")$edges
  selected_node_list <- igraph::as_data_frame(.data, what = "both")$vertices
  all_species <- unique(selected_node_list[, "name"])
  n_edges <- nrow(selected_edge_list)
  
  # loop for groups
  traits_1_edge_list <- selected_edge_list[,1:2]
  for (s in all_species) {
    # get trait for species spp
    trait_1 <- selected_node_list[selected_node_list[, "name"] == s, ][, "group"]
    traits_1_edge_list[selected_edge_list[,1] == as.character(s), 1] <- trait_1
    traits_1_edge_list[selected_edge_list[,2] == as.character(s), 2] <- trait_1
  }
  
  #traits_1_edge_list
  #traits_1_edge_list[traits_1_edge_list[,"to"]=="chiron_l",] # Attention there is a mistake in the trait table of pruning function
  
  # loop for categories
  traits_2_edge_list <- selected_edge_list[,1:2]
  for (s in all_species) {
    # get trait for species spp
    trait_2 <- selected_node_list[selected_node_list[, "name"] == s, ][, "category"]
    traits_2_edge_list[selected_edge_list[,1] == as.character(s), 1] <- trait_2
    traits_2_edge_list[selected_edge_list[,2] == as.character(s), 2] <- trait_2
  }
  
  # loop for variables
  mass_edge_list <- selected_edge_list[,1:2]
  for (s in all_species) {
    # get trait for species spp
    mass.individual <- selected_node_list[selected_node_list[, "name"] == s, ][, "avg.M"]
    mass_edge_list[selected_edge_list[,1] == as.character(s), 1] <- mass.individual
    mass_edge_list[selected_edge_list[,2] == as.character(s), 2] <- mass.individual
  }
  
  # loop for biomasses
  total_edge_list <- selected_edge_list[,1:2]
  for (s in all_species) {
    # get trait for species spp
    mass.individual <- selected_node_list[selected_node_list[, "name"] == s, ][, "avg.M"]
    abundance.individual <- selected_node_list[selected_node_list[, "name"] == s, ][, "avg.N"]
    total_edge_list[selected_edge_list[,1] == as.character(s), 1] <- abundance.individual*mass.individual
    total_edge_list[selected_edge_list[,2] == as.character(s), 2] <- abundance.individual*mass.individual
  }
  
  # merge trait_edge_list with species_edge_list
  merge_edge_traits <- cbind(selected_edge_list[, 1:2],
                             traits_1_edge_list,
                             traits_2_edge_list,
                             mass_edge_list,
                             total_edge_list
  )
  colnames(merge_edge_traits) <- c("from", "to", "trait_from", "trait_to",
                                   "trait_2_from", "trait_2_to",
                                   "mass_consumer", "mass_prey",
                                   "total_mass_predator", "total_mass_prey"
  )
  
  # prepare graph
  prepared_graph <- merge_edge_traits %>%
    graph_from_data_frame()
  
  # The problem is that the selected_node_list does not match the vertex attributes in the graph object
  selected_node_list <- left_join(data.frame(name = as_ids(V(prepared_graph))), selected_node_list, by = "name")
  
  output_graph <- merge_edge_traits %>%
    graph_from_data_frame() %>%
    set_vertex_attr(., c("avg.N"), index = V(.), value = selected_node_list$avg.N) %>%
    set_vertex_attr(., c("avg.M"), index = V(.), value = selected_node_list$avg.M) %>%
    set_vertex_attr(., c("group"), index = V(.), value = selected_node_list$group) %>%
    set_vertex_attr(., c("category"), index = V(.), value = selected_node_list$category)
  
  return(output_graph)
}

# prepare_edge_list_for_pruning for edna datasets which do not have body mass.
prepare_edge_list_for_pruning_edna <- function(.data = NULL)
{
  #.data<- full_metaweb
  if (!is_igraph(.data)) {
    stop("Missing igraph object")
  }
  
  # extract edge list, node list and species list
  selected_edge_list <- igraph::as_data_frame(.data, what = "both")$edges
  selected_node_list <- igraph::as_data_frame(.data, what = "both")$vertices
  all_species <- unique(selected_node_list[, "name"])
  n_edges <- nrow(selected_edge_list)
  
  # loop for groups
  traits_1_edge_list <- selected_edge_list[,1:2]
  for (s in all_species) {
    # get trait for species spp
    trait_1 <- selected_node_list[selected_node_list[, "name"] == s, ][, "group"]
    traits_1_edge_list[selected_edge_list[,1] == as.character(s), 1] <- trait_1
    traits_1_edge_list[selected_edge_list[,2] == as.character(s), 2] <- trait_1
  }
  

  # loop for categories
  traits_2_edge_list <- selected_edge_list[,1:2]
  for (s in all_species) {
    # get trait for species spp
    trait_2 <- selected_node_list[selected_node_list[, "name"] == s, ][, "category"]
    traits_2_edge_list[selected_edge_list[,1] == as.character(s), 1] <- trait_2
    traits_2_edge_list[selected_edge_list[,2] == as.character(s), 2] <- trait_2
  }
  

  
  # merge trait_edge_list with species_edge_list
  merge_edge_traits <- cbind(selected_edge_list[, 1:2],
                             traits_1_edge_list,
                             traits_2_edge_list
  )
  colnames(merge_edge_traits) <- c("from", "to", "trait_from", "trait_to",
                                   "trait_2_from", "trait_2_to"
  )
  
  # prepare graph
  prepared_graph <- merge_edge_traits %>%
    graph_from_data_frame()
  
  # The problem is that the selected_node_list does not match the vertex attributes in the graph object
  selected_node_list <- left_join(data.frame(name = as_ids(V(prepared_graph))), selected_node_list, by = "name")
  
  output_graph <- merge_edge_traits %>%
    graph_from_data_frame() %>%
    set_vertex_attr(., c("avg.N"), index = V(.), value = selected_node_list$avg.N) %>%
    #set_vertex_attr(., c("avg.M"), index = V(.), value = selected_node_list$avg.M) %>%
    set_vertex_attr(., c("group"), index = V(.), value = selected_node_list$group) %>%
    set_vertex_attr(., c("category"), index = V(.), value = selected_node_list$category)
  
  return(output_graph)
}



#' prune_metaweb
#'
#' @param .data
#'
#' @return Pruned metaweb
#' @export
#'
#' @examples prune_full_metaweb(edge_list)
prune_metaweb <- function(.data = NULL, .filter = NULL)
{

  if (!is_igraph(.data)) {
    stop("Missing igraph object")
  }
  
  if (is.null(.filter)) {
    stop("Missing filter selection")
  }
  
  # get edges and nodes from data frame
  input_traits_edge_list <- igraph::as_data_frame(.data, what = "both")$edges
  input_traits_node_list <- igraph::as_data_frame(.data, what = "both")$vertices
  
  if(.filter == 1) {
    # set general interactions rules: first filter
    input_traits_edge_list <- input_traits_edge_list %>%
      subset(., trait_from != "phyto") %>% # remove edges where phyto is a consumer
      subset(., !(trait_from == "zoo" & trait_to == "macro")) %>% # remove edges where zoo preys on macro
      subset(., !(trait_from == "heterotrophic microeukaryot" & trait_to == "zoo")) %>% # remove edges where heterotrophic microeukaryot preys on zoo
      subset(., !(trait_from == "heterotrophic microeukaryot" & trait_to == "macro")) # remove edges where heterotrophic microeukaryot preys on macro
    
  }
  
  if(.filter == 2) {
    # set general interactions rules: first filter
    input_traits_edge_list <- input_traits_edge_list %>%
      subset(., trait_from != "phyto") %>% # remove phyto as consumer
      
      # Microconsumers rules
      subset(., !(trait_from == "zoo" & trait_to == "macro")) %>% # remove zoo  ---> macro (all traits)
      subset(., !(trait_2_from == "zoo_herbivore" & trait_2_to == "zoo_herbivore")) %>% # remove zooplankton grazers ---> grazers
      subset(., !(trait_2_from == "zoo_herbivore" & trait_2_to == "zoo_predator")) %>% # remove zooplankton grazers ---> predator
      subset(., !(trait_2_from == "zoo_herbivore" & trait_2_to == "zoo_omnivore")) %>% # remove zooplankton grazers ---> omnivore
      subset(., !(trait_2_from == "zoo_predator" & trait_to == "phyto")) %>% # remove zooplankton predator ---> phyto
      
      # Macroconsumer rules
      subset(., !(trait_2_from == "macro_herbivore" & trait_2_to == "zoo_herbivore")) %>% # remove macro_herbivore ---> grazers
      subset(., !(trait_2_from == "macro_herbivore" & trait_2_to == "zoo_predator")) %>% # remove macro_herbivore ---> predator
      subset(., !(trait_2_from == "macro_herbivore" & trait_to == "macro")) %>% # remove macro_herbivore ---> macro
      subset(., !(trait_2_from == "macro_predator" & trait_to == "phyto")) %>% # remove macro_predator ---> phyto
      
      # Detritivore and detritus rules
      subset(., (trait_from != "detritus"))  %>% #remove detritus from consumers
      subset(., !(trait_2_from == "macro_detritivore" & trait_to == "macro")) %>%#remove macro from resource of macro_detritivore
      subset(., !(trait_2_from == "macro_detritivore" & trait_to == "zoo")) %>% #remove zoo from resource of macro_detritivore
      subset(., !(trait_2_from == "macro_detritivore" & trait_to == "phyto")) %>%#remove zoo from resource of macro_detritivore
      subset(., !(trait_from == "zoo" & trait_to == "detritus")) %>%#remove zoo from consumer of detritus
      subset(., !(trait_2_from == "macro_herbivore" & trait_to == "detritus")) # remove macro-herbivore from consumer of detritus
    
  }
  
  if(.filter == 3) {
    # set general interactions rules: first filter
    input_traits_edge_list <- input_traits_edge_list %>%
      subset(., trait_from != "phyto") %>% # remove phyto as consumer
      
      # Microconsumers rules
      subset(., !(trait_from == "zoo" & trait_to == "macro")) %>% # remove zoo  ---> macro (all traits)
      
      subset(., !(trait_2_from == "zoo_herbivore" & trait_2_to == "zoo_herbivore")) %>% # remove zooplankton grazers ---> grazers
      subset(., !(trait_2_from == "zoo_herbivore" & trait_2_to == "zoo_predator")) %>% # remove zooplankton grazers ---> predator
      subset(., !(trait_2_from == "zoo_herbivore" & trait_2_to == "zoo_omnivore")) %>% # remove zooplankton grazers ---> omnivore
      
      subset(., !(trait_2_from == "zoo_predator" & trait_to == "phyto")) %>% # remove zooplankton predator ---> phyto
      
      # Macroconsumer rules
      subset(., !(trait_2_from == "macro_herbivore" & trait_2_to == "zoo_herbivore")) %>% # remove macro_herbivore ---> grazers
      subset(., !(trait_2_from == "macro_herbivore" & trait_2_to == "zoo_predator")) %>% # remove macro_herbivore ---> predator
      subset(., !(trait_2_from == "macro_herbivore" & trait_2_to == "zoo_omnivore")) %>% # remove macro_herbivore ---> omnivores
      subset(., !(trait_2_from == "macro_herbivore" & trait_2_to == "macro_herbivore")) %>% # remove macro_herbivore ---> macro_herbivores
      subset(., !(trait_2_from == "macro_herbivore" & trait_2_to == "macro_predator")) %>% # remove macro_herbivore ---> macro_predator
      subset(., !(trait_2_from == "macro_herbivore" & trait_2_to == "macro_omnivore")) %>% # remove macro_herbivore ---> macro_omnivore
      subset(., !(trait_2_from == "macro_herbivore" & trait_2_to == "macro_detritivore")) %>% # remove macro_herbivore ---> macro_detritivore
      subset(., !(trait_2_from == "macro_predator" & trait_to == "phyto")) %>% # remove macro_predator ---> phyto
      
      # Detritivore and detritus rules
      subset(., (trait_from != "detritus")) %>% # remove detritus from consumers
      subset(., !(trait_2_from == "macro_detritivore" & trait_to == "macro")) %>% # remove macro from resource of macro_detritivore
      subset(., !(trait_2_from == "macro_detritivore" & trait_to == "zoo")) %>% # remove zoo from resource of macro_detritivore
      subset(., !(trait_2_from == "macro_detritivore" & trait_to == "phyto")) %>% # remove zoo from resource of macro_detritivore
      subset(., !(trait_from == "zoo" & trait_to == "detritus")) %>% # remove zoo from consumer of detritus
      subset(., !(trait_2_from == "macro_herbivore" & trait_to == "detritus")) # remove macro-herbivore from consumer of detritus
    
  }
  
  # store pruned graph
  pruned_graph <- graph_from_data_frame(input_traits_edge_list, directed = TRUE)
  
  # Check for empty network
  if (length(V(pruned_graph)) == 0) {
    
    output_pruned_graph <- pruned_graph
    
  } else { # graphs with links
    
    # match nodes list with current edge list
    shared_nodes <- intersect(unique(as_ids(V(pruned_graph))), unique(input_traits_node_list$name))
    selected_node_list <- input_traits_node_list[rownames(input_traits_node_list) %in% shared_nodes, ]
    
    # The problem is that the selected_node_list does not match the vertex attributes in the graph object
    selected_node_list <- left_join(data.frame(name = as_ids(V(pruned_graph))), selected_node_list, by = "name")
    
    # reset attributes
    output_pruned_graph <- pruned_graph %>%
      set_vertex_attr(., c("category"), V(.), value = selected_node_list$category) %>%
      set_vertex_attr(., c("group"), V(.), value = selected_node_list$group) %>%
      set_vertex_attr(., c("avg.M"), V(.), value = selected_node_list$avg.M) %>%
      set_vertex_attr(., c("avg.N"), V(.), value = selected_node_list$avg.N)
    
  }
  return(output_pruned_graph)
  
}


#' Prune full metawebs
#'
#' @param .data igraph object
#'
#' @return Pruned metaweb
#' @export
#'
#' @examples prune_full_metaweb(edge_list)
prune_full_metaweb <- function(.data = NULL)
{
  
  if (!is_igraph(.data)) {
    stop("Missing igraph object")
  }
  
  .data <- .data
  # prepare_edge_list_for_pruning(.)
  # get edges and nodes from data frame
  input_traits_edge_list <- igraph::as_data_frame(.data, what = "both")$edges
  input_traits_node_list <- igraph::as_data_frame(.data, what = "both")$vertices
  
  input_traits_edge_list <- input_traits_edge_list %>%
    subset(., trait_from != "phyto") %>% # remove edges where phyto is a consumer
    subset(., !(trait_from == "zoo" & trait_to == "macro")) %>% # remove edges where zoo preys on macro
    subset(., !(trait_from == "macro" & trait_to == "phyto")) # remove edges where zoo preys on macro
  
  # store pruned graph
  pruned_graph <- graph_from_data_frame(input_traits_edge_list, directed = TRUE)
  
  # match nodes list with current list
  shared_nodes <- intersect(unique(as_ids(V(pruned_graph))), unique(input_traits_node_list$name))
  selected_node_list <- input_traits_node_list[rownames(input_traits_node_list) %in% shared_nodes, ]
  
  # reset attributes
  output_pruned_graph <- pruned_graph %>%
    set_vertex_attr(., c("avg.N"), index = V(.), value = selected_node_list$avg.N) %>%
    set_vertex_attr(., c("avg.M"), index = V(.), value = selected_node_list$avg.M) %>%
    set_vertex_attr(., c("group"), index = V(.), value = selected_node_list$group) %>%
    set_vertex_attr(., c("category"), index = V(.), value = selected_node_list$category)
  
  return(output_pruned_graph)
}