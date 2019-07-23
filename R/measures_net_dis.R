#' Netdis between all graph pairs using provided Centred Graphlet Counts
#' @param centred_graphlet_counts List containing Centred Graphlet Counts for 
#' all graphs being compared
#' @param graphlet_size The size of graphlets to use for the Netdis calculation
#' (only counts for graphlets of the specified size will be used). The size of
#' a graphlet is the number of nodes it contains.
#' @return Pairwise Netdis statistics between graphs calculated using centred 
#' counts for graphlets of the specified size
#' @export
netdis_for_all_graphs <- function(
  centred_graphlet_counts, graphlet_size, mc.cores = getOption("mc.cores", 2L)) {
  comp_spec <- cross_comparison_spec(centred_graphlet_counts)
  # NOTE: mcapply only works on unix-like systems with system level forking 
  # capability. This means it will work on Linux and OSX, but not Windows.
  # For now, we just revert to single threaded operation on Windows
  # TODO: Look into using the parLappy function on Windows
  if(.Platform$OS.type != "unix") {
    # Force cores to 1 if system is not unix-like as it will not support 
    # forking
    mc.cores = 1
  }
  netdis <- purrr::simplify(parallel::mcmapply(function(index_a, index_b) {netdis(
    centred_graphlet_counts[[index_a]], centred_graphlet_counts[[index_b]], 
    graphlet_size = graphlet_size)
  }, comp_spec$index_a, comp_spec$index_b, SIMPLIFY = FALSE))
  list(netdis = netdis, comp_spec = comp_spec)
}

#' Netdis
#' 
#' Calculate Netdis statistic between two graphs from their Centred Graphlet
#' Counts (generated using \code{netdis_centred_graphlet_counts}).
#' @param centred_graphlet_counts1 Centred Graphlet Counts for graph 1
#' @param centred_graphlet_counts2 Centred Graphlet Counts for graph 2
#' @param graphlet_size The size of graphlets to use for the Netdis calculation
#' (only counts for graphlets of the specified size will be used). The size of
#' a graphlet is the number of nodes it contains.
#' @return Netdis statistic calculated using centred counts for graphlets of 
#' the specified size
#' @export
netdis <- function(centred_graphlet_counts1, centred_graphlet_counts2, 
                   graphlet_size)
{
  # Select subset of centred counts corresponding to graphlets of the 
  # specified size
  ids <- graphlet_ids_for_size(graphlet_size)
  counts1 <- centred_graphlet_counts1[ids]
  counts2 <- centred_graphlet_counts2[ids]
  
  # Calculate normalising constant
  norm_const <- sum(counts1^2 / sqrt(counts1^2 + counts2^2),na.rm = TRUE) *
    sum(counts2^2 / sqrt(counts1^2 + counts2^2),na.rm = TRUE)
  # Calculate intermediate "netD" statistic that falls within range -1..1
  netds2 <- (1/sqrt(norm_const)) * sum((counts1 * counts2) / sqrt(counts1^2 + counts2^2),na.rm = TRUE)
  # Calculate corresponding "netd" Netdis statistic that falls within range 0..1
  0.5 * (1 - netds2)
} 

#' Scaled graphlet count for ego-networks
#' 
#' Calculates graphlet counts for the n-step ego-network of each node in a graph, 
#' scaled by dividing the graphlet counts for each ego-network by the total
#' number of possible groupings of nodes in the ego-network with the same number
#' of nodes as each graphlet. This scaling factor is choose(n, k), where n is the
#' number of nodes in the ego-network and k is the number of nodes in the graphlet.
#' @param graph A connected, undirected, simple graph as an \code{igraph} object. 
#' @param max_graphlet_size Determines the maximum size of graphlets to count. 
#' Only graphlets containing up to \code{max_graphlet_size} nodes will be counted.
#' @param neighbourhood_size The number of steps from the source node to include
#' nodes for each ego-network.
#' @param return_ego_networks If \code{TRUE}, return ego-networks alongside 
#' graphlet counts to enable further processing. 
#' @return If \code{return_ego_networks = FALSE}, returns an RxC matrix 
#' containing counts of each graphlet (columns, C) for each ego-network in the 
#' input graph (rows, R). Columns are labelled with graphlet IDs and rows are 
#' labelled with the ID of the central node in each ego-network (if nodes in the
#' input graph are labelled). If \code{return_ego_networks = TRUE}, returns a
#' list with the following elements:
#' \itemize{
#'   \item \code{graphlet_counts}: A matrix containing graphlet counts for each 
#'   ego-network in the input graph as described above.
#'   \item \code{ego_networks}: The ego-networks of the query graph.
#' }
#' @export
count_graphlets_ego_scaled <- function(
  graph, max_graphlet_size, neighbourhood_size, return_ego_networks = FALSE) {
  # Calculate ego-network graphlet counts, also returning the ego networks for
  # use later in function
  ego_data <- 
    count_graphlets_ego(graph, max_graphlet_size = max_graphlet_size, 
                        neighbourhood_size = neighbourhood_size, 
                        return_ego_networks = TRUE)
  ego_graphlet_counts <- ego_data$graphlet_counts
  ego_networks <- ego_data$ego_networks
  # Scale ego-network graphlet counts by dividing by total number of k-tuples in
  # ego-network (where k is graphlet size)
  ego_graphlet_tuples <- 
    count_graphlet_tuples_ego(ego_networks, max_graphlet_size = max_graphlet_size)
  ego_graphlet_counts <- scale_graphlet_count(ego_graphlet_counts, ego_graphlet_tuples)
  # Return either graphlet counts, or graphlet counts and ego_networks
  if(return_ego_networks) {
    return(list(graphlet_counts = ego_graphlet_counts, 
                ego_networks = ego_networks))
  } else {
    return(ego_graphlet_counts)
  }
}

#' Generate Netdis centred graphlets counts by subtracting expected counts
#' 
#' @param graph A connected, undirected, simple graph as an \code{igraph} object.
#' @param max_graphlet_size Determines the maximum size of graphlets to count. 
#' Only graphlets containing up to \code{max_graphlet_size} nodes will be counted.
#' @param neighbourhood_size The number of steps from the source node to include
#' nodes for each ego-network. 
#' @param expected_ego_count_fn A function for generating expected ego-network
#' graphlet counts for a graph. This function should take a connected, 
#' undirected, simple graph as an \code{igraph} object for its only argument. 
#' Where \code{expected_ego_count_fn} is specific to particular values of 
#' \code{max_graphlet_size} or \code{neighbourhood_size}, care should be taken 
#' to ensure that the values of these parameters passed to this function are
#' consistent with those used when creating \code{expected_ego_count_fn}.
#' @return A vector with centred counts for each graphlet type
#' @export 
netdis_centred_graphlet_counts <- function(
  graph, max_graphlet_size, neighbourhood_size, expected_ego_count_fn = 
    NULL) {
  # Get centred counts for each ego network
  centred_counts <- netdis_centred_graphlet_counts_ego(
    graph, max_graphlet_size, neighbourhood_size, expected_ego_count_fn)
  # Sum centred counts over ego-networks
  apply(centred_counts, MARGIN = 2, FUN = sum)
}


#' TODO: Remove @export prior to publishing
#' @export
netdis_centred_graphlet_counts_ego <- function(
  graph, max_graphlet_size, neighbourhood_size, expected_ego_count_fn = NULL,
  min_ego_nodes = 3, min_ego_edges = 1) {
    # Get unscaled ego-network graphlet counts
    res <- count_graphlets_ego(
      graph, max_graphlet_size = max_graphlet_size,
      neighbourhood_size = neighbourhood_size, return_ego_networks = TRUE)
    actual_counts = res$graphlet_counts
    ego_networks <- res$ego_networks
    
    # Drop ego-networks that don't have the minimum number of nodes or edges
    drop_index <- purrr::simplify(purrr::map(ego_networks, function(g) { 
      (igraph::vcount(g) < min_ego_nodes) | (igraph::ecount(g) < min_ego_edges)
    }))
    actual_counts <- actual_counts[!drop_index,]
    ego_networks <- ego_networks[!drop_index]
    
    # Centre these counts by subtracting the expected counts
    if(is.null(expected_ego_count_fn)) {
      centred_counts = actual_counts
    } else {
      centred_counts <- actual_counts - expected_ego_count_fn(graph)
    }
    centred_counts
}

#' Generate Netdis expected graphlet count function
#' 
#' Generates a function to calculate expected ego-network graphlet counts for 
#' query graphs based on the statistics of a provided reference graph.
#' 
#' Generates graphlet counts for all ego-networks in the supplied reference graph 
#' and then averages these graphlet counts over density bins to generate
#' density-dependent reference graphlet counts. Prior to averaging, the graphlet
#' counts are scaled in a size-dependent manner to permit ego-networks with 
#' similar densities but different sizes to be averaged together.
#' 
#' Returns a function that uses the density-dependent reference graphlet 
#' counts to generate expected graphlet counts for all ego-networks in a query
#' network. When doing so, it matches ego-networks to reference counts by
#' density and reverses the scaling that was applied to the original reference
#' counts in order to allow averaging across ego-networks with similar density
#' but different numbers of nodes.
#' @param graph A connected, undirected, simple reference graph as an 
#' \code{igraph} object. 
#' @param max_graphlet_size Determines the maximum size of graphlets to count. 
#' Only graphlets containing up to \code{max_graphlet_size} nodes will be counted.
#' @param neighbourhood_size The number of steps from the source node to include
#' node in ego-network.
#' @return A function taking a connected, undirected, simple query graph as an 
#' \code{igraph} object and returning an RxC matrix containing the expected 
#' counts of each graphlet (columns, C) for each ego-network in the query graph 
#' (rows, R). Columns are labelled with graphlet IDs and rows are labelled with 
#' the ID of the central node in each ego-network (if nodes in the query graph 
#' are labelled)
#' @export
netdis_expected_graphlet_counts_ego_fn <- function(
  graph, max_graphlet_size, neighbourhood_size,
  min_ego_nodes = 3, min_ego_edges = 1, 
  min_bin_count = 5, num_bins = 100) {
  
  # Calculate the scaled graphlet counts for all ego networks in the reference
  # graph, also returning the ego networks themselves in order to calculate
  # their densities
  res <- count_graphlets_ego_scaled(
    graph, max_graphlet_size, neighbourhood_size, return_ego_networks = TRUE)
  scaled_graphlet_counts = res$graphlet_counts
  ego_networks <- res$ego_networks
  
  # Drop ego-networks that don't have the minimum number of nodes or edges
  # JACK  - why not put this in make_named_ego_graph? i.e. when generating ego networks in first place
  drop_index <- purrr::simplify(purrr::map(ego_networks, function(g) { 
    (igraph::vcount(g) < min_ego_nodes) | (igraph::ecount(g) < min_ego_edges)
  }))
  scaled_graphlet_counts <- scaled_graphlet_counts[!drop_index,]
  ego_networks <- ego_networks[!drop_index]
  
  # Get ego-network densities
  densities <- purrr::simplify(purrr::map_dbl(ego_networks, igraph::edge_density))
  
  # Adaptively bin ego-network densities
  binned_densities <- binned_densities_adaptive(
    densities, min_counts_per_interval = min_bin_count, num_intervals = num_bins)
  
  # Average graphlet counts across density bins
  density_binned_graphlet_counts <- mean_density_binned_graphlet_counts(
    scaled_graphlet_counts, binned_densities$interval_indexes)
  
  # Return a partially applied function with the key reference graph information
  # built-in
  purrr::partial(
    netdis_expected_graphlet_counts_ego,
    max_graphlet_size = max_graphlet_size,
    neighbourhood_size = neighbourhood_size,
    min_ego_nodes = min_ego_nodes, 
    min_ego_edges = min_ego_edges,
    density_breaks = binned_densities$breaks,
    density_binned_reference_counts = density_binned_graphlet_counts)
}

#' INTERNAL FUNCTION - Do not call directly
#' 
#' Used by \code{netdis_expected_graphlet_counts_ego_fn} to 
#' generate a function for calculating expected ego-network graphlet counts
#' from the statistics of a provided reference graph.
#' Temporarily accessible during development. 
#' TODO: Remove @export prior to publishing
#' @export
netdis_expected_graphlet_counts_ego <- function(
  graph, max_graphlet_size, neighbourhood_size,
  density_breaks, density_binned_reference_counts,
  min_ego_nodes = 3, min_ego_edges = 1) {
  # Generate ego-networks for query graph
  ego_networks <- make_named_ego_graph(graph, neighbourhood_size)
  # Drop ego-networks that don't have the minimum number of nodes or edges
  drop_index <- purrr::simplify(purrr::map(ego_networks, function(g) { 
    (igraph::vcount(g) < min_ego_nodes) | (igraph::ecount(g) < min_ego_edges)
  }))
  ego_networks <- ego_networks[!drop_index]
  # Map over query graph ego-networks, using reference graph statistics to 
  # calculate expected graphlet counts for each ego-network.
  expected_graphlet_counts <- 
    purrr::map(ego_networks, netdis_expected_graphlet_counts,
               max_graphlet_size = max_graphlet_size,
               density_breaks = density_breaks,
               density_binned_reference_counts = density_binned_reference_counts)
  names(expected_graphlet_counts) <- names(ego_networks)
  # Simplify list to array
  t(simplify2array(expected_graphlet_counts))
}

#' INTERNAL FUNCTION - Do not call directly
#' 
#' Used by \code{netdis_expected_graphlet_counts_ego} to 
#' calculate expected graphlet counts for a query graph ego-network from the 
#' statistics of a provided reference graph.
#' Temporarily accessible during development. 
#' TODO: Remove @export prior to publishing
#' @export
netdis_expected_graphlet_counts <- function(
  graph, max_graphlet_size, density_breaks, density_binned_reference_counts) {
  # Look up average scaled graphlet counts for graphs of similar density
  # in the reference graph
  query_density <- igraph::edge_density(graph)
  matched_density_index <- interval_index(query_density, density_breaks)
  matched_reference_counts <- density_binned_reference_counts[matched_density_index,]
  # Scale reference counts by multiplying the reference count for each graphlet
  # by the number of possible sets of k nodes in the query graph, where k is the
  # number of nodes in the graphlet
  matched_reference_counts * count_graphlet_tuples(graph, max_graphlet_size)
}

#' INTERNAL FUNCTION - Do not call directly
#' 
#' Used by \code{netdis_expected_graphlet_counts_ego_fn} to 
#' generate a function for calculating expected ego-network graphlet counts
#' from the statistics of a provided reference graph.
#' Temporarily accessible during development. 
#' TODO: Remove @export prior to publishing
#' @export
mean_density_binned_graphlet_counts <- function(
  graphlet_counts, density_interval_indexes) {
  # The ego network graphlet counts are an E x G matrix with rows (E) representing
  # ego networks and columns (G) representing graphlets. We want to calculate
  # the mean count for each graphlet / density bin combination, so we will
  # use tapply to average counts for each graphlet across density bins, using
  # apply to map this operation over graphlets
  mean_density_binned_graphlet_counts <- 
    apply(graphlet_counts, MARGIN = 2, function(gc) {
      tapply(gc, INDEX = density_interval_indexes, FUN = mean)})
}


#' @export
zeros_to_ones <- function(v) {
  zero_index <- which(v == 0)
  v[zero_index] <- 1
  v
}

#' @export
scale_graphlet_count <- function(graphlet_count, graphlet_tuples) {
  # Avoid divide by zero errors by replacing all zeros with ones in the
  # divisor
  graphlet_count / zeros_to_ones(graphlet_tuples)
}

#' @export
count_graphlet_tuples_ego <- function(ego_networks, max_graphlet_size) {
  graphlet_tuple_counts <- 
    t(simplify2array(purrr::map(ego_networks, count_graphlet_tuples, 
             max_graphlet_size = max_graphlet_size)))
  graphlet_tuple_counts
}

#' @export
count_graphlet_tuples <- function(graph, max_graphlet_size) {
  graph_node_count <- igraph::vcount(graph)
  graphlet_key <- graphlet_key(max_graphlet_size)
  graphlet_node_counts <- graphlet_key$node_count
  graphlet_tuple_counts <- choose(graph_node_count, graphlet_node_counts)
  graphlet_tuple_counts <- stats::setNames(graphlet_tuple_counts, graphlet_key$id)
  graphlet_tuple_counts
}



