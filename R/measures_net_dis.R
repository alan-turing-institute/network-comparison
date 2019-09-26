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
#' @param binning_fn Function used to bin ego network densities. Takes densities
#' as its single argument, and returns a named list including keys \code{breaks}
#' (list of bin edges) and \code{interval_indexes} (density bin index for each
#' ego network). (Default: \code{binned_densities_adaptive} with
#' \code{min_counts_per_interval = 5} and \code{num_intervals = 100}).
#' @param bin_counts_fn Function used to calculate expected graphlet counts in
#' each density bin. Takes \code{graphlet_counts}, \code{interval_indexes}
#' (bin indexes), \code{ego_networks} and \code{max_graphlet_size} as arguments.
#' (Default: \code{density_binned_counts} with \code{agg_fn = mean} and
#' \code{scale_fn = scale_graphlet_counts_ego}, which mirrors the
#' approach used in the original netdis paper).
#' @param exp_counts_fn Function used to map from binned reference counts to
#' expected counts for each graphlet in each ego network of the query graphs.
#' Takes \code{ego_networks}, \code{density_bin_breaks},
#' \code{binned_graphlet_counts}, and \code{max_graphlet_size} as arguments.
#' (Default: \code{netdis_expected_graphlet_counts_per_ego} with
#' \code{scale_fn = count_graphlet_tuples}, which mirrors the approach used in
#' the original netdis paper).
#' @return Netdis statistics between graph_1 and graph_2 for graphlet sizes
#' up to and including max_graphlet_size
#' @export
netdis_one_to_one <- function(graph_1, graph_2,
                              ref_graph,
                              max_graphlet_size = 4,
                              neighbourhood_size = 2,
                              min_ego_nodes = 3,
                              min_ego_edges = 1,
                              binning_fn = purrr::partial(
                                binned_densities_adaptive,
                                min_counts_per_interval = 5,
                                num_intervals = 100),
                              bin_counts_fn = purrr::partial(
                                density_binned_counts,
                                agg_fn = mean,
                                scale_fn = scale_graphlet_counts_ego),
                              exp_counts_fn = purrr::partial(
                                netdis_expected_graphlet_counts_per_ego,
                                scale_fn = count_graphlet_tuples)) {

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
#' @param binning_fn Function used to bin ego network densities. Takes densities
#' as its single argument, and returns a named list including keys \code{breaks}
#' (list of bin edges) and \code{interval_indexes} (density bin index for each
#' ego network). (Default: \code{binned_densities_adaptive} with
#' \code{min_counts_per_interval = 5} and \code{num_intervals = 100}).
#' @param bin_counts_fn Function used to calculate expected graphlet counts in
#' each density bin. Takes \code{graphlet_counts}, \code{interval_indexes}
#' (bin indexes), \code{ego_networks} and \code{max_graphlet_size} as arguments.
#' (Default: \code{density_binned_counts} with \code{agg_fn = mean} and
#' \code{scale_fn = scale_graphlet_counts_ego}, which mirrors the
#' approach used in the original netdis paper).
#' @param exp_counts_fn Function used to map from binned reference counts to
#' expected counts for each graphlet in each ego network of the query graphs.
#' Takes \code{ego_networks}, \code{density_bin_breaks},
#' \code{binned_graphlet_counts}, and \code{max_graphlet_size} as arguments.
#' (Default: \code{netdis_expected_graphlet_counts_per_ego} with
#' \code{scale_fn = count_graphlet_tuples}, which mirrors the approach used in
#' the original netdis paper).
#' @return Netdis statistics between graph_1 and graph_2 for graphlet sizes
#' up to and including max_graphlet_size
#' @export
netdis_one_to_many <- function(graph_1, graphs_compare,
                              ref_graph,
                              max_graphlet_size = 4,
                              neighbourhood_size = 2,
                              min_ego_nodes = 3,
                              min_ego_edges = 1,
                              binning_fn = purrr::partial(
                                binned_densities_adaptive,
                                min_counts_per_interval = 5,
                                num_intervals = 100),
                              bin_counts_fn = purrr::partial(
                                density_binned_counts,
                                agg_fn = mean,
                                scale_fn = scale_graphlet_counts_ego),
                              exp_counts_fn = purrr::partial(
                                netdis_expected_graphlet_counts_per_ego,
                                scale_fn = count_graphlet_tuples)) {

  # bundle graph_1 and graphs_compare to one vector, with
  # graph_1 at start as needed for netdis_many_to_many call
  graphs <- append(graphs_compare, list(graph_1 = graph_1), after = 0)

  # calculate netdis
  result <- netdis_many_to_many(
    graphs,
    ref_graph,
    comparisons = "one-to-many",
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
#' @param binning_fn Function used to bin ego network densities. Takes densities
#' as its single argument, and returns a named list including keys \code{breaks}
#' (list of bin edges) and \code{interval_indexes} (density bin index for each
#' ego network). (Default: \code{binned_densities_adaptive} with
#' \code{min_counts_per_interval = 5} and \code{num_intervals = 100}).
#' @param bin_counts_fn Function used to calculate expected graphlet counts in
#' each density bin. Takes \code{graphlet_counts}, \code{interval_indexes}
#' (bin indexes), \code{ego_networks} and \code{max_graphlet_size} as arguments.
#' (Default: \code{density_binned_counts} with \code{agg_fn = mean} and
#' \code{scale_fn = scale_graphlet_counts_ego}, which mirrors the
#' approach used in the original netdis paper).
#' @param exp_counts_fn Function used to map from binned reference counts to
#' expected counts for each graphlet in each ego network of the query graphs.
#' Takes \code{ego_networks}, \code{density_bin_breaks},
#' \code{binned_graphlet_counts}, and \code{max_graphlet_size} as arguments.
#' (Default: \code{netdis_expected_graphlet_counts_per_ego} with
#' \code{scale_fn = count_graphlet_tuples}, which mirrors the approach used in
#' the original netdis paper).
#' @return Netdis statistics between query graphs for graphlet sizes
#' up to and including max_graphlet_size.
#' @export
netdis_many_to_many <- function(graphs,
                                ref_graph,
                                comparisons = "many-to-many",
                                max_graphlet_size = 4,
                                neighbourhood_size = 2,
                                min_ego_nodes = 3,
                                min_ego_edges = 1,
                                binning_fn = purrr::partial(
                                  binned_densities_adaptive,
                                  min_counts_per_interval = 5,
                                  num_intervals = 100),
                                bin_counts_fn = purrr::partial(
                                  density_binned_counts,
                                  agg_fn = mean,
                                  scale_fn = scale_graphlet_counts_ego),
                                exp_counts_fn = purrr::partial(
                                  netdis_expected_graphlet_counts_per_ego,
                                  scale_fn = count_graphlet_tuples)) {
  ## ------------------------------------------------------------------------
  # Get ego networks for query graphs
  ego_networks <- purrr::map(
    graphs, make_named_ego_graph,
    order = neighbourhood_size,
    min_ego_nodes = min_ego_nodes,
    min_ego_edges = min_ego_edges
  )

  rm(graphs)
  
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
    
    rm(ref_graph)

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

    # Calculate expected graphlet counts (using ref
    # graph ego network density bins)
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
    densities <- purrr::map(graphlet_counts,
                            ego_network_density)

    # bin ref ego-network densities
    binned_densities <- purrr::map(densities,
                                   binning_fn)

    # extract bin breaks and indexes from binning results
    ego_density_bin_breaks <- purrr::map(binned_densities,
                                          function(x) {
                                            x$breaks
                                          })
    ego_density_bin_indexes <- purrr::map(binned_densities,
                                          function(x) {
                                            x$interval_indexes
                                          })


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

#' INTERNAL FUNCTION - Do not call directly
#'
#' Used by \code{netdis_expected_graphlet_counts_ego_fn} to
#' generate a function for calculating expected ego-network graphlet counts
#' from the statistics of a provided reference graph.
#'
#' @param ego_networks The number of steps from the source node to include
#' node in ego-network.
#' @param density_breaks Density values defining bin edges.
#' @param density_binned_reference_counts Reference network graphlet counts for
#' each density bin.
#' @param max_graphlet_size Determines the maximum size of graphlets to count.
#' Only graphlets containing up to \code{max_graphlet_size} nodes are counted.
#' @param scale_fn Optional function to scale calculated expected counts, taking
#' \code{graph} and \code{max_graphlet_size} as arguments, and returning a scale
#' factor that the looked up \code{density_binned_reference_counts} values will
#' be multiplied by.
#'
#' #' Temporarily accessible during development.
#' TODO: Remove @export prior to publishing
#' @export
netdis_expected_graphlet_counts_per_ego <- function(
                                              ego_networks,
                                              density_breaks,
                                              density_binned_reference_counts,
                                              max_graphlet_size,
                                              scale_fn=NULL) {


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
  expected_graphlet_counts <- t(simplify2array(expected_graphlet_counts))
  colnames(expected_graphlet_counts) <- graphlet_key(max_graphlet_size)$id
  expected_graphlet_counts
}

#' INTERNAL FUNCTION - Do not call directly
#'
#' Used by \code{netdis_expected_graphlet_counts_ego} to
#' calculate expected graphlet counts for a query graph
#' ego-network from the statistics of a provided reference
#' graph.
#' @param graph A connected, undirected, simple reference graph as an
#' \code{igraph} object.
#' @param max_graphlet_size Determines the maximum size of graphlets to count.
#' Only graphlets containing up to \code{max_graphlet_size} nodes are counted.
#' @param density_breaks Density values defining bin edges.
#' @param density_binned_reference_counts Reference network graphlet counts for
#' each density bin.
#' @param scale_fn Optional function to scale calculated expected counts, taking
#' \code{graph} and \code{max_graphlet_size} as arguments, and returning a scale
#' factor that the looked up \code{density_binned_reference_counts} values will
#' be multiplied by.
#' Temporarily accessible during development.
#' TODO: Remove @export prior to publishing
#' @export
netdis_expected_graphlet_counts <- function(graph,
                                            max_graphlet_size,
                                            density_breaks,
                                            density_binned_reference_counts,
                                            scale_fn=NULL) {

  # Look up average scaled graphlet counts for graphs of similar density
  # in the reference graph
  query_density <- igraph::edge_density(graph)
  matched_density_index <- interval_index(query_density, density_breaks)

  matched_reference_counts <-
    density_binned_reference_counts[matched_density_index, ]
  
  if (!is.null(scale_fn)) {
    # Scale reference counts e.g. by multiplying the
    # reference count for each graphlet by the number
    # of possible sets of k nodes in the query graph,
    # where k is the number of nodes in the graphlet.
    matched_reference_counts <- matched_reference_counts *
                                scale_fn(graph, max_graphlet_size)
  }
  matched_reference_counts
}

#' INTERNAL FUNCTION - Do not call directly
#'
#' Used by \code{netdis_expected_graphlet_counts_ego_fn} to
#' generate a function for calculating expected ego-network graphlet counts
#' from the statistics of a provided reference graph.
#' @param graphlet_counts Graphlet counts for a number of ego_networks.
#' @param density_interval_indexes Density bin index for
#' each ego network.
#' @param agg_fn Function to aggregate counts in each bin
#' (default \code{agg_fn = mean}).
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
    dim(mean_density_binned_graphlet_counts) <-
      c(1, length(mean_density_binned_graphlet_counts))

    colnames(mean_density_binned_graphlet_counts) <-
      colnames(graphlet_counts)
  }

  mean_density_binned_graphlet_counts
}

#' For case where don't want to use binning, return a single bin which covers
#' the full range of possible density values.
#' @param densities Ego network density values (only used to return
#' a list of indexes of the required length.)
#' @export
single_density_bin <- function(densities) {
  list(densities = densities,
       interval_indexes = rep(1, length(densities)),
       breaks = c(0, 1))
}

#' INTERNAL FUNCTION - Do not call directly
#'
#' Used to calculate expected graphlet counts for each density bin.
#' @param graphlet_counts Graphlet counts for a number of ego_networks.
#' @param density_interval_indexes Density bin index for
#' each ego network.
#' @param agg_fn Function to aggregate counts in each bin
#' (default \code{agg_fn = mean}).
#' @param scale_fn Optional function to apply a transformation
#' to graphlet_counts, must have arguments graphlet_counts,
#' ego_networks and max_graphlet_size.
#' @param ego_networks Optionally passed and used by scale_fn.
#' @param max_graphlet_size Optionally passed and used by scale_fn.
#' @export
density_binned_counts <- function(graphlet_counts,
                                  density_interval_indexes,
                                  agg_fn = mean,
                                  scale_fn = NULL,
                                  ego_networks = NULL,
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


#' Calculate expected counts in density bins using
#' geometric poisson (Polya-Aeppli) approximation.
#' @param graphlet_counts Graphlet counts for a number of ego_networks.
#' @param density_interval_indexes Density bin index for
#' each ego network.
#' @param max_graphlet_size Determines the maximum size of graphlets
#' included in graphlet_counts.
#' @export
density_binned_counts_gp <- function(graphlet_counts,
                                     density_interval_indexes,
                                     max_graphlet_size) {

  mean_binned_graphlet_counts <- mean_density_binned_graphlet_counts(
    graphlet_counts,
    density_interval_indexes)

  exp_counts_bin <- function(bin_idx) {
    counts <- graphlet_counts[density_interval_indexes == bin_idx, ]
    means <- mean_binned_graphlet_counts[bin_idx, ]

    mean_sub_counts <- sweep(counts, 2, means)

    Vd_sq <- colSums(mean_sub_counts^2) / (nrow(mean_sub_counts) - 1)
    theta_d <- 2 * means / (Vd_sq + means)

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
  expected_counts_bin[is.nan(expected_counts_bin)] <- 0

  expected_counts_bin
}


#' Replace zero values in a vector with ones. Used by
#' \code{scale_graphlet_count} to prevent divide by
#' zero errors.
#' @param v A vector.
#' TODO remove export
#' @export
zeros_to_ones <- function(v) {
  zero_index <- which(v == 0)
  v[zero_index] <- 1
  v
}


#' Divide graphlet counts by pre-computed scaling factor from
#' \code{count_graphlet_tuples} output.
#' @param graphlet_count Pre-computed graphlet counts.
#' @param graphlet_tuples Pre-computed \code{count_graphlet_tuples} output.
#' @export
scale_graphlet_count <- function(graphlet_count, graphlet_tuples) {
  # Avoid divide by zero errors by replacing all zeros with ones in the
  # divisor
  graphlet_count / zeros_to_ones(graphlet_tuples)
}


#' Run count_graphlet_tuples across pre-computed ego networks.
#' @param ego_networks Pre-generated ego networks for an input graph.
#' @param max_graphlet_size Determines the maximum size of graphlets included
#' in the tuple counts.
#' @export
count_graphlet_tuples_ego <- function(graphlet_counts, max_graphlet_size) {
  graphlet_tuple_counts <-
    t(apply(graphlet_counts, 1,
            count_graphlet_tuples, max_graphlet_size = max_graphlet_size))
  
  graphlet_tuple_counts
}


#' Calculate edge density for a single graph.
#' @param graphlet_counts Vector of pre-calculated graphlet, edge and node 
#' counts. Must have named items "N" (node counts) and "G0" (edge counts).
#' @export
density_from_counts <- function(graphlet_counts) {
  graphlet_counts["G0"] / choose(graphlet_counts["N"], 2)
}

#' Calculate ego network edge densities.
#' @param graphlet_counts Matrix of pre-generated graphlet, edge and node counts
#' (columns) for each ego network (rows). Columns must include "N" (node counts)
#' and "G0" (edge counts).
#' @export
ego_network_density <- function(graphlet_counts) {
  apply(graphlet_counts, 1, density_from_counts)
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
scale_graphlet_counts_ego <- function(graphlet_counts,
                                      max_graphlet_size) {
  ego_graphlet_tuples <- count_graphlet_tuples_ego(
    graphlet_counts,
    max_graphlet_size = max_graphlet_size
  )

  scaled_graphlet_counts <- scale_graphlet_count(
    graphlet_counts,
    ego_graphlet_tuples
  )

  return(scaled_graphlet_counts)
}


#' For each graphlet calculate the number of possible sets of k nodes in the
#' query graph, where k is the number of nodes in the graphlet.
#'
#' @param graph A connected, undirected, simple graph as an \code{igraph}
#' object.
#' @param max_graphlet_size Determines the maximum size of graphlets included
#' in the tuple counts.
#' @export
count_graphlet_tuples <- function(graph_graphlet_counts, max_graphlet_size) {
  # extract node counts from graph_graphlet_counts
  N <- graph_graphlet_counts["N"]
  
  graphlet_key <- graphlet_key(max_graphlet_size)
  graphlet_node_counts <- graphlet_key$node_count
  
  graphlet_tuple_counts <- choose(N, graphlet_node_counts)
  
  graphlet_tuple_counts <- stats::setNames(
    graphlet_tuple_counts,
    graphlet_key$id
  )
  
  # add node counts back to object
  graphlet_tuple_counts <- c(N, graphlet_tuple_counts)
}
