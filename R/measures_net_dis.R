#' Netdis between two graphs
#' @param graph_1 First query graph
#' @param graph_2 Second query graph
#' @param ref_graph Reference graph
#' @param max_graphlet_size Generate graphlets up to this size
#' @param neighbourhood_size Ego network neighbourhood size
#' @param min_ego_nodes Filter ego networks which have fewer
#' than min_ego_nodes nodes
#' @param min_ego_edges Filter ego networks which have fewer
#' than min_ego_edges edges
#' @param min_bin_count Minimum number of ego networks in each density bin
#' @param num_bins Number of density bins to generate
#' @return Netdis statistics between graph_1 and graph_2 for graphlet sizes
#' up to and including max_graphlet_size
#' @export
netdis_one_to_one <- function(graph_1, graph_2,
                              ref_graph,
                              max_graphlet_size = 4,
                              neighbourhood_size = 2,
                              min_ego_nodes = 3,
                              min_ego_edges = 1,
                              binning_fn = purrr::partial(binned_densities_adaptive, min_counts_per_interval = 5, num_intervals = 100),
                              bin_counts_fn = purrr::partial(density_binned_counts, agg_fn = mean, scale_fn = scale_graphlet_counts_ego),
                              exp_counts_fn = purrr::partial(netdis_expected_graphlet_counts_per_ego, scale_fn=count_graphlet_tuples)) {
  
  # bundle graphs into one vector with format needed for
  # netdis many-to-many
  graphs <- list(graph_1 = graph_1, graph_2 = graph_2)
  
  # calculate netdis
  result <- netdis_many_to_many(
    graphs,
    ref_graph,
    max_graphlet_size = 4,
    neighbourhood_size = 2,
    min_ego_nodes = 3,
    min_ego_edges = 1,
    binning_fn = binning_fn,
    exp_counts_fn = exp_counts_fn
  )
  
  # extract netdis statistics from list returned by netdis_many_to_many
  result$netdis[, 1]
}

#' Netdis comparisons between one graph and many other graphs
#' @param graph_1 query graph - this graph will be compared with
#' all graphs in graphs_compare
#' @param graphs_compare graphs graph_1 will be compared with
#' @param ref_graph Reference graph
#' @param max_graphlet_size Generate graphlets up to this size
#' @param neighbourhood_size Ego network neighbourhood size
#' @param min_ego_nodes Filter ego networks which have fewer
#' than min_ego_nodes nodes
#' @param min_ego_edges Filter ego networks which have fewer
#' than min_ego_edges edges
#' @param min_bin_count Minimum number of ego networks in each density bin
#' @param num_bins Number of density bins to generate
#' @return Netdis statistics between graph_1 and graph_2 for graphlet sizes
#' up to and including max_graphlet_size
#' @export
netdis_one_to_many <- function(graph_1, graphs_compare,
                              ref_graph,
                              max_graphlet_size = 4,
                              neighbourhood_size = 2,
                              min_ego_nodes = 3,
                              min_ego_edges = 1,
                              binning_fn = purrr::partial(binned_densities_adaptive, min_counts_per_interval = 5, num_intervals = 100),
                              bin_counts_fn = purrr::partial(density_binned_counts, agg_fn = mean, scale_fn = scale_graphlet_counts_ego),
                              exp_counts_fn = purrr::partial(netdis_expected_graphlet_counts_per_ego, scale_fn=count_graphlet_tuples)) {
  
  # bundle graph_1 and graphs_compare to one vector, with
  # graph_1 at start as needed for netdis_many_to_many call
  graphs <- append(graphs_compare, list(graph_1=graph_1), after=0)

  # calculate netdis
  result <- netdis_many_to_many(
    graphs,
    ref_graph,
    comparisons = 'one-to-many',
    max_graphlet_size = 4,
    neighbourhood_size = 2,
    min_ego_nodes = 3,
    min_ego_edges = 1,
    binning_fn = binning_fn,
    bin_counts_fn = bin_counts_fn,
    exp_counts_fn = exp_counts_fn
  )
  
  # restructure netdis_many_to_many output
  colnames(result$netdis) <- result$comp_spec$name_b
  result$netdis
}


#' Netdis between all graph pairs
#' @param graphs Query graphs
#' @param ref_graph Reference graph
#' @param comparisons Which comparisons to perform between graphs.
#' Can be "many-to-many" (all pairwise combinations) or "one-to-many"
#' (compare first graph in graphs to all other graphs.)
#' @param max_graphlet_size Generate graphlets up to this size
#' @param neighbourhood_size Ego network neighbourhood size
#' @param min_ego_nodes Filter ego networks which have fewer
#' than min_ego_nodes nodes
#' @param min_ego_edges Filter ego networks which have fewer
#' than min_ego_edges edges
#' @param min_bin_count Minimum number of ego networks in each density bin
#' @param num_bins Number of density bins to generate
#' @return Netdis statistics between graph_1 and graph_2 for graphlet sizes
#' up to and including max_graphlet_size
#' @export
netdis_many_to_many <- function(graphs,
                                ref_graph,
                                comparisons = 'many-to-many',
                                max_graphlet_size = 4,
                                neighbourhood_size = 2,
                                min_ego_nodes = 3,
                                min_ego_edges = 1,
                                binning_fn = purrr::partial(binned_densities_adaptive,
                                                            min_counts_per_interval = 5,
                                                            num_intervals = 100),
                                bin_counts_fn = purrr::partial(density_binned_counts,
                                                               agg_fn = mean,
                                                               scale_fn = scale_graphlet_counts_ego),
                                exp_counts_fn = purrr::partial(netdis_expected_graphlet_counts_per_ego,
                                                               scale_fn=count_graphlet_tuples)) {
  ## ------------------------------------------------------------------------
  # Get ego networks for query graphs
  ego_networks <- purrr::map(
    graphs, make_named_ego_graph,
    order = neighbourhood_size, 
    min_ego_nodes = min_ego_nodes, 
    min_ego_edges = min_ego_edges
  )
  
  ## ------------------------------------------------------------------------
  # Count graphlets for ego networks in query graphs
  graphlet_counts <- purrr::map(
    ego_networks,
    ego_to_graphlet_counts,
    max_graphlet_size = max_graphlet_size
  )
  
  ## ------------------------------------------------------------------------
  # Case where expected counts calculated using a reference network
  if (!is.null(ref_graph)) {
    # Get ego networks
    ego_ref <- make_named_ego_graph(
      ref_graph, 
      order = neighbourhood_size, 
      min_ego_nodes = min_ego_nodes, 
      min_ego_edges = min_ego_edges
    )
    
    # Get ego network graphlet counts
    graphlet_counts_ref <- ego_to_graphlet_counts(
      ego_ref,
      max_graphlet_size = max_graphlet_size
    )
    
    # Get ego-network densities
    densities_ref <- ego_network_density(ego_ref)
    
    # bin ref ego-network densities
    binned_densities <- binning_fn(densities_ref)
    
    ref_ego_density_bins <- binned_densities$breaks
  
    # Average ref graphlet counts across density bins
    ref_binned_graphlet_counts <- bin_counts_fn(
                                    graphlet_counts_ref, 
                                    binned_densities$interval_indexes,
                                    ego_networks = ego_ref,
                                    max_graphlet_size = max_graphlet_size
                                  )
    
    # Calculate expected graphlet counts (using ref graph ego network density bins)
    exp_graphlet_counts <- purrr::map(
      ego_networks,
      exp_counts_fn,
      density_breaks = ref_ego_density_bins,
      density_binned_reference_counts = ref_binned_graphlet_counts,
      max_graphlet_size = max_graphlet_size
    )
    
  ## ------------------------------------------------------------------------
  } else {
    # Case where expected counts calculated using query networks
    
    # Get ego-network densities
    densities <- purrr::map(ego_networks,
                            ego_network_density)
    
    # bin ref ego-network densities
    binned_densities <- purrr::map(densities,
                                   binning_fn)
    
    # extract bin breaks and indexes from binning results
    ego_density_bin_breaks <- purrr::map(binned_densities,
                                          function(x) {x$breaks})
    ego_density_bin_indexes <- purrr::map(binned_densities,
                                          function(x) {x$interval_indexes})
    
    
    # Calculate expected counts in each bin
    binned_graphlet_counts <- mapply(bin_counts_fn,
                                     graphlet_counts,
                                     ego_density_bin_indexes,
                                     max_graphlet_size = max_graphlet_size,
                                     SIMPLIFY = FALSE)
    
    # Calculate expected graphlet counts for each ego network
    exp_graphlet_counts <- mapply(exp_counts_fn,
                                  ego_networks,
                                  ego_density_bin_breaks,
                                  binned_graphlet_counts,
                                  max_graphlet_size = max_graphlet_size,
                                  SIMPLIFY = FALSE)
  }
  
  ## ------------------------------------------------------------------------  
  # Centre graphlet counts by subtracting expected counts
  centred_graphlet_counts <- mapply("-", graphlet_counts, exp_graphlet_counts)
  
  ## ------------------------------------------------------------------------
  # Sum centred graphlet counts across all ego networks
  sum_graphlet_counts <- lapply(centred_graphlet_counts, colSums)
  
  ## ------------------------------------------------------------------------
  # Generate pairwise comparisons
  comp_spec <- cross_comparison_spec(sum_graphlet_counts, how = comparisons)
  
  ## ------------------------------------------------------------------------
  # Calculate netdis statistics
  results <- parallel::mcmapply(
      function(index_a, index_b) {
        netdis_uptok(
          sum_graphlet_counts[[index_a]], 
          sum_graphlet_counts[[index_b]],
          max_graphlet_size = max_graphlet_size
        )
      },
      comp_spec$index_a,
      comp_spec$index_b,
      SIMPLIFY = TRUE)
  
  
  list(netdis = results, comp_spec = comp_spec)
  
}

#' Netdis between all graph pairs using provided Centred Graphlet Counts
#' @param centred_graphlet_counts List containing Centred Graphlet Counts for
#' all graphs being compared
#' @param graphlet_size The size of graphlets to use for the Netdis calculation
#' (only counts for graphlets of the specified size will be used). The size of
#' a graphlet is the number of nodes it contains.
#' @return Pairwise Netdis statistics between graphs calculated using centred
#' counts for graphlets of the specified size
#' @export
netdis_for_all_graphs <- function(centred_graphlet_counts,
                                  graphlet_size,
                                  mc.cores = getOption("mc.cores", 2L)) {
  comp_spec <- cross_comparison_spec(centred_graphlet_counts)
  # NOTE: mcapply only works on unix-like systems with system level forking
  # capability. This means it will work on Linux and OSX, but not Windows.
  # For now, we just revert to single threaded operation on Windows
  # TODO: Look into using the parLappy function on Windows
  if (.Platform$OS.type != "unix") {
    # Force cores to 1 if system is not unix-like as it will not support
    # forking
    mc.cores <- 1
  }
  netdis <- purrr::simplify(parallel::mcmapply(function(index_a, index_b) {
    netdis(
      centred_graphlet_counts[[index_a]], centred_graphlet_counts[[index_b]],
      graphlet_size = graphlet_size
    )
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
                   graphlet_size) {
  # Select subset of centred counts corresponding to graphlets of the
  # specified size
  ids <- graphlet_ids_for_size(graphlet_size)
  counts1 <- centred_graphlet_counts1[ids]
  counts2 <- centred_graphlet_counts2[ids]

  # Calculate normalising constant
  norm_const <- sum(counts1^2 / sqrt(counts1^2 + counts2^2), na.rm = TRUE) *
    sum(counts2^2 / sqrt(counts1^2 + counts2^2), na.rm = TRUE)
  # Calculate intermediate "netD" statistic that falls within range -1..1
  netds2 <- (1 / sqrt(norm_const)) *
    sum((counts1 * counts2) /
      sqrt(counts1^2 + counts2^2), na.rm = TRUE)
  # Calculate corresponding "netd" Netdis statistic that falls within range 0..1
  0.5 * (1 - netds2)
}

#' Netdis - graphlets up to max_graphlet_size
#'
#' Calculate Netdis statistic between two graphs from their Centred Graphlet
#' Counts (generated using \code{netdis_centred_graphlet_counts}).
#' @param centred_graphlet_counts1 Centred Graphlet Counts for graph 1
#' @param centred_graphlet_counts2 Centred Graphlet Counts for graph 2
#' @param max_graphlet_size max graphlet size to calculate Netdis for.
#' The size of a graphlet is the number of nodes it contains. Netdis is
#' calculated for all graphlets from size 3 to size max_graphlet_size.
#' @return Netdis statistic calculated using centred counts for graphlets of
#' the specified size
#' @export
netdis_uptok <- function(centred_graphlet_counts1, centred_graphlet_counts2,
                         max_graphlet_size) {
  if ((max_graphlet_size > 5) | (max_graphlet_size < 3)) {
    stop("max_graphlet_size must be 3, 4 or 5.")
  }

  netdis_statistics <- purrr::map(3:max_graphlet_size,
    netdis,
    centred_graphlet_counts1 = centred_graphlet_counts1,
    centred_graphlet_counts2 = centred_graphlet_counts2
  )

  netdis_statistics <- simplify2array(netdis_statistics)

  names(netdis_statistics) <-
    sapply(
      "netdis",
      paste,
      3:max_graphlet_size,
      sep = ""
    )

  netdis_statistics
}

#' Scaled graphlet count for ego-networks
#'
#' Calculates graphlet counts for the n-step ego-network of each node in
#' a graph, scaled by dividing the graphlet counts for each ego-network by the
#' total number of possible groupings of nodes in the ego-network with the same
#' number of nodes as each graphlet. This scaling factor is choose(n, k),
#' where n is the number of nodes in the ego-network and k is the number of
#' nodes in the graphlet.
#' @param graph A connected, undirected, simple graph as an \code{igraph}
#'  object.
#' @param max_graphlet_size Determines the maximum size of graphlets to count.
#' Only graphlets containing up to \code{max_graphlet_size} nodes are counted.
#' @param neighbourhood_size The number of steps from the source node to include
#' nodes for each ego-network.
#' @param min_ego_nodes Only ego networks with at least \code{min_ego_nodes}
#' nodes are returned.
#' @param min_ego_edges Only ego networks with at least \code{min_ego_edges}
#' edges are returned.
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
count_graphlets_ego_scaled <- function(graph,
                                       max_graphlet_size,
                                       neighbourhood_size,
                                       min_ego_nodes = 3,
                                       min_ego_edges = 1,
                                       return_ego_networks = FALSE) {

  # Calculate ego-network graphlet counts, also returning the ego networks for
  # use later in function
  ego_data <-
    count_graphlets_ego(graph,
      max_graphlet_size = max_graphlet_size,
      neighbourhood_size = neighbourhood_size,
      min_ego_nodes = min_ego_nodes,
      min_ego_edges = min_ego_edges,
      return_ego_networks = TRUE
    )
  ego_graphlet_counts <- ego_data$graphlet_counts
  ego_networks <- ego_data$ego_networks

  # Scale ego-network graphlet counts by dividing by total number of k-tuples in
  # ego-network (where k is graphlet size)
  ego_graphlet_tuples <- count_graphlet_tuples_ego(
    ego_networks,
    max_graphlet_size = max_graphlet_size
  )
  ego_graphlet_counts <- scale_graphlet_count(
    ego_graphlet_counts,
    ego_graphlet_tuples
  )

  # Return either graphlet counts, or graphlet counts and ego_networks
  if (return_ego_networks) {
    return(list(
      graphlet_counts = ego_graphlet_counts,
      ego_networks = ego_networks
    ))
  } else {
    return(ego_graphlet_counts)
  }
}

#' Generate Netdis centred graphlets counts by subtracting expected counts
#'
#' @param graph A connected, undirected, simple graph as an
#' \code{igraph} object.
#' @param max_graphlet_size Determines the maximum size of graphlets to count.
#' Only graphlets containing up to \code{max_graphlet_size} nodes
#' will be counted.
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
netdis_centred_graphlet_counts <- function(graph,
                                           max_graphlet_size,
                                           neighbourhood_size,
                                           expected_ego_count_fn = NULL) {
  # Get centred counts for each ego network
  centred_counts <- netdis_centred_graphlet_counts_ego(
    graph,
    max_graphlet_size,
    neighbourhood_size,
    expected_ego_count_fn
  )
  # Sum centred counts over ego-networks
  apply(centred_counts, MARGIN = 2, FUN = sum)
}


#' TODO: Remove @export prior to publishing
#' @export
netdis_centred_graphlet_counts_ego <- function(graph,
                                               max_graphlet_size,
                                               neighbourhood_size,
                                               expected_ego_count_fn = NULL,
                                               min_ego_nodes = 3,
                                               min_ego_edges = 1) {
  # Get unscaled ego-network graphlet counts
  res <- count_graphlets_ego(
    graph,
    max_graphlet_size = max_graphlet_size,
    neighbourhood_size = neighbourhood_size,
    min_ego_nodes = min_ego_nodes,
    min_ego_edges = min_ego_edges,
    return_ego_networks = TRUE
  )

  actual_counts <- res$graphlet_counts

  # Centre these counts by subtracting the expected counts
  if (is.null(expected_ego_count_fn)) {
    centred_counts <- actual_counts
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
#' Generates graphlet counts for all ego-networks in the supplied
#' reference graph and then averages these graphlet counts over density bins to
#' generate density-dependent reference graphlet counts. Prior to averaging,
#' the graphlet counts are scaled in a size-dependent manner to permit
#' ego-networks with similar densities but different sizes to be averaged
#' together.
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
#' Only graphlets containing up to \code{max_graphlet_size} nodes are counted.
#' @param neighbourhood_size The number of steps from the source node to include
#' node in ego-network.
#' @return A function taking a connected, undirected, simple query graph as an
#' \code{igraph} object and returning an RxC matrix containing the expected
#' counts of each graphlet (columns, C) for each ego-network in the query graph
#' (rows, R). Columns are labelled with graphlet IDs and rows are labelled with
#' the ID of the central node in each ego-network (if nodes in the query graph
#' are labelled)
#' @export
netdis_expected_graphlet_counts_ego_fn <- function(graph,
                                                   max_graphlet_size,
                                                   neighbourhood_size,
                                                   min_ego_nodes = 3,
                                                   min_ego_edges = 1,
                                                   min_bin_count = 5,
                                                   num_bins = 100,
                                                   scale_fn = NULL) {

  # Calculate the scaled graphlet counts for all ego networks in the reference
  # graph, also returning the ego networks themselves in order to calculate
  # their densities
  res <- count_graphlets_ego_scaled(
    graph,
    max_graphlet_size,
    neighbourhood_size,
    min_ego_nodes = min_ego_nodes,
    min_ego_edges = min_ego_edges,
    return_ego_networks = TRUE
  )

  scaled_graphlet_counts <- res$graphlet_counts
  ego_networks <- res$ego_networks

  # Get ego-network densities
  densities <- purrr::simplify(
    purrr::map_dbl(ego_networks, igraph::edge_density)
  )

  # Adaptively bin ego-network densities
  binned_densities <- binned_densities_adaptive(
    densities,
    min_counts_per_interval = min_bin_count,
    num_intervals = num_bins
  )

  # Average graphlet counts across density bins
  density_binned_graphlet_counts <- mean_density_binned_graphlet_counts(
    scaled_graphlet_counts,
    binned_densities$interval_indexes
  )

  # Return a partially applied function with the key reference graph information
  # built-in
  purrr::partial(
    netdis_expected_graphlet_counts_ego,
    max_graphlet_size = max_graphlet_size,
    neighbourhood_size = neighbourhood_size,
    min_ego_nodes = min_ego_nodes,
    min_ego_edges = min_ego_edges,
    density_breaks = binned_densities$breaks,
    density_binned_reference_counts = density_binned_graphlet_counts,
    scale_fn = scale_fn
  )
}

#' INTERNAL FUNCTION - Do not call directly
#'
#' Used by \code{netdis_expected_graphlet_counts_ego_fn} to
#' generate a function for calculating expected ego-network graphlet counts
#' from the statistics of a provided reference graph.
#' Temporarily accessible during development.
#' TODO: Remove @export prior to publishing
#' @export
netdis_expected_graphlet_counts_ego <- function(graph,
                                                max_graphlet_size,
                                                neighbourhood_size,
                                                density_breaks,
                                                density_binned_reference_counts,
                                                min_ego_nodes = 3,
                                                min_ego_edges = 1,
                                                scale_fn=NULL) {
  
  #print("netdis_expected_graphlet_counts_ego")
  #print(density_binned_reference_counts)
  
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
      density_binned_reference_counts = density_binned_reference_counts,
      scale_fn=scale_fn
    )
  names(expected_graphlet_counts) <- names(ego_networks)
  # Simplify list to array
  t(simplify2array(expected_graphlet_counts))
}

#' INTERNAL FUNCTION - Do not call directly
#'
#' JACK To follow through logic of paper steps, wanted to pass
#' ego networks to the function, not the input query graph
#' (as in netdis_expected_graphlet_counts_ego above).
#'
#' Used by \code{netdis_expected_graphlet_counts_ego_fn} to
#' generate a function for calculating expected ego-network graphlet counts
#' from the statistics of a provided reference graph.
#' Temporarily accessible during development.
#' TODO: Remove @export prior to publishing
#' @export
netdis_expected_graphlet_counts_per_ego <- function(ego_networks,
                                                    density_breaks,
                                                    density_binned_reference_counts,
                                                    max_graphlet_size,
                                                    scale_fn=NULL) {
  
  
  #print("netdis_expected_graphlet_counts_per_ego")
  #print(density_binned_reference_counts)
  
  # Map over query graph ego-networks, using reference graph statistics to
  # calculate expected graphlet counts for each ego-network.
  expected_graphlet_counts <-
    purrr::map(ego_networks, netdis_expected_graphlet_counts,
      max_graphlet_size = max_graphlet_size,
      density_breaks = density_breaks,
      density_binned_reference_counts = density_binned_reference_counts,
      scale_fn = scale_fn
    )
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
netdis_expected_graphlet_counts <- function(graph,
                                            max_graphlet_size,
                                            density_breaks,
                                            density_binned_reference_counts,
                                            scale_fn=NULL) {
  
  #print("netdis_expected_graphlet_counts")
  #print(density_binned_reference_counts)
  
  # Look up average scaled graphlet counts for graphs of similar density
  # in the reference graph
  query_density <- igraph::edge_density(graph)
  matched_density_index <- interval_index(query_density, density_breaks)

  matched_reference_counts <-
    density_binned_reference_counts[matched_density_index, ]
  
  if (!is.null(scale_fn)) {
    # Scale reference counts e.g. by multiplying the reference count for each graphlet
    # by the number of possible sets of k nodes in the query graph, where k is the
    # number of nodes in the graphlet
    matched_reference_counts <- matched_reference_counts * scale_fn(graph, max_graphlet_size)
  }
  
  matched_reference_counts
}

#' INTERNAL FUNCTION - Do not call directly
#'
#' Used by \code{netdis_expected_graphlet_counts_ego_fn} to
#' generate a function for calculating expected ego-network graphlet counts
#' from the statistics of a provided reference graph.
#' Temporarily accessible during development.
#' TODO: Remove @export prior to publishing
#' @export
mean_density_binned_graphlet_counts <- function(graphlet_counts,
                                                density_interval_indexes,
                                                agg_fn = mean) {
  # The ego network graphlet counts are an E x G matrix with rows (E)
  # representing ego networks and columns (G) representing graphlets. We want
  # to calculate the mean count for each graphlet / density bin combination,
  # so we will use tapply to average counts for each graphlet across density
  # bins, using apply to map this operation over graphlets
  mean_density_binned_graphlet_counts <-
    apply(graphlet_counts, MARGIN = 2, function(gc) {
      tapply(gc, INDEX = density_interval_indexes, FUN = agg_fn)
    })
  
  # if only 1 bin (i.e. no binning) will be left with a 1D list.
  # convert it into a 2D list.
  if (is.null(dim(mean_density_binned_graphlet_counts))) {
    dim(mean_density_binned_graphlet_counts) <- c(1, length(mean_density_binned_graphlet_counts))
    colnames(mean_density_binned_graphlet_counts) <- colnames(graphlet_counts)
  }
  
  mean_density_binned_graphlet_counts
}

#' For case where don't want to use binning, return
#' a single bin which covers full range of possible
#' densities.
#' @export
single_density_bin <- function(densities) {

  binned_densities <- list(densities = densities,
                           interval_indexes = rep(1, length(densities)),
                           breaks = c(0, 1))
}

#' INTERNAL FUNCTION - Do not call directly
#'
#' Used to
#' generate a function for calculating expected graphlet counts in each
#' density bin.
#' @param agg_fn Function to aggregate counts in each bin (default \code{agg_fn = mean}).
#' @param scale_fn Optional function to apply a transformation to graphlet_counts, must
#' have arguments graphlet_counts, ego_networks and max_graphlet_size.
#' @param ego_networks Optionally passed and used by scale_fn.
#' @param max_graphlet_size Optionally passed and used by scale_fn.
#' Temporarily accessible during development.
#' TODO: Remove @export prior to publishing
#' @export
density_binned_counts <- function(graphlet_counts, density_interval_indexes,
                                  agg_fn = mean,
                                  scale_fn = NULL, ego_networks = NULL,
                                  max_graphlet_size = NULL) {
  
  if (!is.null(scale_fn)) {
    # Scale ego-network graphlet counts e.g.
    # by dividing by total number of k-tuples in 
    # ego-network (where k is graphlet size)
    graphlet_counts <- scale_fn(graphlet_counts,
                                ego_networks = ego_networks,
                                max_graphlet_size = max_graphlet_size)
  }
  
  mean_density_binned_graphlet_counts(graphlet_counts,
                                      density_interval_indexes,
                                      agg_fn = agg_fn)
  
}


#' Calculate expected counts in density bins using geometric poisson (Polya-Aeppli) approximation
#' @export
density_binned_counts_gp <- function(graphlet_counts,
                                     density_interval_indexes,
                                     max_graphlet_size) {
  
  mean_binned_graphlet_counts <- mean_density_binned_graphlet_counts(
    graphlet_counts, 
    density_interval_indexes)
  
  exp_counts_bin <- function(bin_idx) {
    counts <- graphlet_counts[density_interval_indexes == bin_idx, ]
    means <- mean_binned_graphlet_counts[bin_idx,]
    
    mean_sub_counts <- sweep(counts, 2, means)
    
    Vd_sq <- colSums(mean_sub_counts^2)/(nrow(mean_sub_counts)-1)
    theta_d <- 2*means / (Vd_sq + means)
    
    exp_counts_dk <- vector()
    for (k in 2:max_graphlet_size) {
      graphlet_idx <- graphlet_ids_for_size(k)
      
      lambda_dk <- (1 / length(graphlet_idx)) * 
        sum(
          2 * means[graphlet_idx]^2 /
            (Vd_sq[graphlet_idx] + means[graphlet_idx])
        )
      
      exp_counts_dk <- append(exp_counts_dk,
                              lambda_dk / theta_d[graphlet_idx])
    }
    
    exp_counts_dk
  }
  
  nbins <- length(unique(density_interval_indexes))
  expected_counts_bin <- t(mapply(exp_counts_bin, bin_idx = 1:nbins))
  
  # deal with NAs caused by bins with zero counts for a graphlet
  expected_counts_bin[is.nan(expected_counts_bin)] = 0
  
  expected_counts_bin
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
      max_graphlet_size = max_graphlet_size
    )))
  graphlet_tuple_counts
}

#' @export
ego_network_density <- function(ego_networks) {
  densities <- purrr::simplify(purrr::map_dbl(
    ego_networks,
    igraph::edge_density
  ))

  return(densities)
}



#' Scale graphlet counts for an ego network by the n choose k possible
#' choices of k nodes in that ego-network, where n is the number of nodes
#' in the ego network and k is the number of nodes in the graphlet.
#'
#' @param ego_networks Pre-generated ego networks for an input graph.
#' @param graphlet_counts Pre-calculated graphlet counts for each ego_network.
#' @param max_graphlet_size Determines the maximum size of graphlets included
#' in graphlet_counts.
#' @return scaled graphlet counts.
#' @export
scale_graphlet_counts_ego <- function(ego_networks, graphlet_counts,
                                      max_graphlet_size) {
  ego_graphlet_tuples <- count_graphlet_tuples_ego(
    ego_networks,
    max_graphlet_size = max_graphlet_size
  )

  scaled_graphlet_counts <- scale_graphlet_count(
    graphlet_counts,
    ego_graphlet_tuples
  )

  return(scaled_graphlet_counts)
}


#' @export
count_graphlet_tuples <- function(graph, max_graphlet_size) {
  graph_node_count <- igraph::vcount(graph)
  graphlet_key <- graphlet_key(max_graphlet_size)
  graphlet_node_counts <- graphlet_key$node_count
  graphlet_tuple_counts <- choose(graph_node_count, graphlet_node_counts)
  graphlet_tuple_counts <- stats::setNames(
    graphlet_tuple_counts,
    graphlet_key$id
  )
  graphlet_tuple_counts
}
