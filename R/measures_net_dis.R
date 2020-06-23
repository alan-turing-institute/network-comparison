#' Netdis between two graphs
#'
#' Calculates the different variants of the network dissimilarity statistic Netdis between two graphs. The variants currently supported are Netdis using a gold-standard network, Netdis using no expecations (\code{ref_graph = 0}), and Netdis using a Geometric Poisson  approximation for the expectation (\code{ref_graph = NULL}).
#' 
#' 
#' @param graph_1 A simplified igraph graph object.
#'
#' @param graph_2 A simplified igraph graph object.
#' 
#' @param graphlet_counts_1 Pre-generated graphlet counts for the first query
#' graph. If the \code{graphlet_counts_1} argument is defined then
#' \code{graph_1} will not be used.
#'
#' @param graphlet_counts_2 Pre-generated graphlet counts for the second query
#' graph. If the \code{graphlet_counts_2} argument is defined then
#' \code{graph_2} will not be used.
#'
#' @param ref_graph Controls how expected counts are calculated. Either:
#' 1) A numeric value - used as a constant expected counts value for all query
#' graphs (DEFAULT: 0).
#' 2) A simplified \code{igraph} object - used as a reference graph from which
#' expected counts are calculated for all query graphs.
#' 3) NULL - Expected counts will be calculated based on the properties of the
#' query graphs themselves.
#' 
#' @param graphlet_counts_ref Pre-generated reference graphlet counts. If the
#' \code{graphlet_counts_ref} argument is defined then \code{ref_graph} will not
#' be used.
#'
#' @param max_graphlet_size Generate graphlets up to this size.
#'
#' @param neighbourhood_size Ego network neighborhood size.
#'
#' @param min_ego_nodes Filter ego networks which have fewer
#' than min_ego_nodes nodes.
#'
#' @param min_ego_edges Filter ego networks which have fewer
#' than min_ego_edges edges.
#'
#' @param binning_fn Function used to bin ego network densities. Takes densities
#' as its single argument, and returns a named list including keys \code{breaks}
#' (list of bin edges) and \code{interval_indexes} (density bin index for each
#' ego network). (Default: \code{binned_densities_adaptive} with
#' \code{min_counts_per_interval = 5} and \code{num_intervals = 100}).
#'
#' @param bin_counts_fn Function used to calculate expected graphlet counts in
#' each density bin. Takes \code{graphlet_counts}, \code{interval_indexes}
#' (bin indexes) and \code{max_graphlet_size} as arguments.
#' (Default: \code{density_binned_counts} with \code{agg_fn = mean} and
#' \code{scale_fn = scale_graphlet_counts_ego}, which mirrors the
#' approach used in the original netdis paper).
#'
#' @param exp_counts_fn Function used to map from binned reference counts to
#' expected counts for each graphlet in each ego network of the query graphs.
#' Takes \code{ego_networks}, \code{density_bin_breaks},
#' \code{binned_graphlet_counts}, and \code{max_graphlet_size} as arguments.
#' (Default: \code{netdis_expected_counts} with
#' \code{scale_fn = count_graphlet_tuples}, which mirrors the approach used in
#' the original netdis paper).
#'
#' @return Netdis statistics between graph_1 and graph_2 for graphlet sizes
#' up to and including max_graphlet_size.
#'
#' @examples
#' require(netdist)
#' require(igraph)
#' #Set source directory for Virus PPI graph edge files stored in the netdist package.
#' source_dir <- system.file(file.path("extdata", "VRPINS"), package = "netdist")
#' # Load query graphs as igraph objects
#' graph_1 <- read_simple_graph(file.path(source_dir, "EBV.txt"),format = "ncol")
#' graph_2 <- read_simple_graph(file.path(source_dir, "ECL.txt"),format = "ncol")
#' 
#' #Netdis variant using the Geometric Poisson approximation to remove the background expectation of each network.
#' netdis_one_to_one(graph_1= graph_1, graph_2= graph_2,  ref_graph = NULL) #This option will focus on detecting more general and global discrepancies between the ego-network structures.
#' 
#' #Comparing the networks via their observed ego counts without centering them (equivalent to using expectation equal to zero). This option, will focus on detecting small discrepancies.
#' netdis_one_to_one(graph_1= graph_1, graph_2= graph_2,  ref_graph = 0)
#' 
#' # Example of the use of netdis with a reference graph.This option will focus on detecting discrepancies between the networks relative to the ego-network structure of the reference network / gold-standard.
#' # Two lattice networks of different sizes are used for this example. 
#'  goldstd_1 <- graph.lattice(c(8,8)) #A reference net
#'  goldstd_2 <- graph.lattice(c(44,44)) #A reference net
#'  
#'  netdis_one_to_one(graph_1= graph_1, graph_2= graph_2,  ref_graph = goldstd_1)
#'  netdis_one_to_one(graph_1= graph_1, graph_2= graph_2,  ref_graph = goldstd_2)
#'  
#'  
#'  #Imputing pre-calculated subgraph counts instead of subgraphs.
#'  
#'  props_1 <- count_graphlets_ego(graph = graph_1)
#'  props_2 <- count_graphlets_ego(graph = graph_2)
#'  props_goldstd_1 <- count_graphlets_ego(graph = goldstd_1)
#'  props_goldstd_2 <- count_graphlets_ego(graph = goldstd_2)
#'  
#' #Netdis Geometric-Poisson.
#' netdis_one_to_one(graphlet_counts_1= props_1,graphlet_counts_2= props_2, ref_graph = NULL)
#' 
#' #Netdis Zero Expectation.
#' netdis_one_to_one(graphlet_counts_1= props_1,graphlet_counts_2= props_2, ref_graph = 0)
#' 
#' #Netdis using gold-standard network
#' netdis_one_to_one(graphlet_counts_1= props_1,graphlet_counts_2= props_2, graphlet_counts_ref = props_goldstd_1)
#' netdis_one_to_one(graphlet_counts_1= props_1,graphlet_counts_2= props_2, graphlet_counts_ref = props_goldstd_2)
#' @export
netdis_one_to_one <- function(graph_1 = NULL,
                              graph_2 = NULL,
                              ref_graph = 0,
                              max_graphlet_size = 4,
                              neighbourhood_size = 2,
                              min_ego_nodes = 3,
                              min_ego_edges = 1,
                              binning_fn = purrr::partial(
                                binned_densities_adaptive,
                                min_counts_per_interval = 5,
                                num_intervals = 100
                              ),
                              bin_counts_fn = purrr::partial(
                                density_binned_counts,
                                agg_fn = mean,
                                scale_fn = scale_graphlet_counts_ego
                              ),
                              exp_counts_fn = purrr::partial(
                                netdis_expected_counts,
                                scale_fn = count_graphlet_tuples
                              ),
                              graphlet_counts_1 = NULL,
                              graphlet_counts_2 = NULL,
                              graphlet_counts_ref= NULL) {
  
  ## ------------------------------------------------------------------------
  # Check arguments
  if (is.null(graph_1) & is.null(graphlet_counts_1)) {
    stop("One of graph_1 and graphlet_counts_1 must be supplied.")
  }
  if (is.null(graph_2) & is.null(graphlet_counts_2)) {
    stop("One of graph_2 and graphlet_counts_2 must be supplied.")
  }
  ## ------------------------------------------------------------------------
  # Generate graphlet counts and bundle them into named list with format needed
  # for netdis_many_to_many.
  
  if (is.null(graphlet_counts_1)) {
    graphlet_counts_1 <- count_graphlets_ego(
      graph_1,
      max_graphlet_size = max_graphlet_size,
      neighbourhood_size = neighbourhood_size,
      min_ego_nodes = min_ego_nodes,
      min_ego_edges = min_ego_edges,
      return_ego_networks = FALSE
    )
  }
  rm(graph_1)
  
  if (is.null(graphlet_counts_2)) {
    graphlet_counts_2 <- count_graphlets_ego(
      graph_2,
      max_graphlet_size = max_graphlet_size,
      neighbourhood_size = neighbourhood_size,
      min_ego_nodes = min_ego_nodes,
      min_ego_edges = min_ego_edges,
      return_ego_networks = FALSE
    )
  }
  rm(graph_2)
  
  graphlet_counts <- list(
    graph_1 = graphlet_counts_1,
    graph_2 = graphlet_counts_2
  )
  
  if(!is.null(ref_graph)){
    if (!is.numeric(ref_graph) & is.null(graphlet_counts_ref)) {
      graphlet_counts_ref <- count_graphlets_ego(
        ref_graph,
        max_graphlet_size = max_graphlet_size,
        neighbourhood_size = neighbourhood_size,
        min_ego_nodes = min_ego_nodes,
        min_ego_edges = min_ego_edges,
        return_ego_networks = FALSE
      )
      ref_graph <- NULL
    }
  }
  ## ------------------------------------------------------------------------
  # calculate netdis
  result <- netdis_many_to_many(
    graphs = NULL,
    ref_graph = ref_graph,
    max_graphlet_size = 4,
    neighbourhood_size = 2,
    min_ego_nodes = 3,
    min_ego_edges = 1,
    binning_fn = binning_fn,
    bin_counts_fn = bin_counts_fn,
    exp_counts_fn = exp_counts_fn,
    graphlet_counts = graphlet_counts,
    graphlet_counts_ref = graphlet_counts_ref
  )
  
  ## ------------------------------------------------------------------------
  # extract netdis statistics from list returned by netdis_many_to_many
  result$netdis[, 1]
}





#' Netdis comparisons between one graph and many other graphs.
#'
#' @param graph_1 Query graph - this graph will be compared with
#' all graphs in graphs_compare. A simplified igraph graph object.
#'
#' @param graphs_compare Graphs graph_1 will be compared with. A named list of
#' simplified igraph graph objects.
#'
#' @param ref_graph Controls how expected counts are calculated. Either:
#' 1) A numeric value - used as a constant expected counts value for all query
#' graphs (DEFAULT: 0).
#' 2) A simplified \code{igraph} object - used as a reference graph from which
#' expected counts are calculated for all query graphs.
#' 3) NULL - Expected counts will be calculated based on the properties of the
#' query graphs themselves.
#'
#' @param max_graphlet_size Generate graphlets up to this size.
#'
#' @param neighbourhood_size Ego network neighbourhood size.
#'
#' @param min_ego_nodes Filter ego networks which have fewer
#' than min_ego_nodes nodes.
#'
#' @param min_ego_edges Filter ego networks which have fewer
#' than min_ego_edges edges.
#'
#' @param binning_fn Function used to bin ego network densities. Takes densities
#' as its single argument, and returns a named list including keys \code{breaks}
#' (list of bin edges) and \code{interval_indexes} (density bin index for each
#' ego network). (Default: \code{binned_densities_adaptive} with
#' \code{min_counts_per_interval = 5} and \code{num_intervals = 100}).
#'
#' @param bin_counts_fn Function used to calculate expected graphlet counts in
#' each density bin. Takes \code{graphlet_counts}, \code{interval_indexes}
#' (bin indexes) and \code{max_graphlet_size} as arguments.
#' (Default: \code{density_binned_counts} with \code{agg_fn = mean} and
#' \code{scale_fn = scale_graphlet_counts_ego}, which mirrors the
#' approach used in the original netdis paper).
#'
#' @param exp_counts_fn Function used to map from binned reference counts to
#' expected counts for each graphlet in each ego network of the query graphs.
#' Takes \code{ego_networks}, \code{density_bin_breaks},
#' \code{binned_graphlet_counts}, and \code{max_graphlet_size} as arguments.
#' (Default: \code{netdis_expected_counts} with
#' \code{scale_fn = count_graphlet_tuples}, which mirrors the approach used in
#' the original netdis paper).
#'
#' @param graphlet_counts_1 Pre-generated graphlet counts for the first query
#' graph. If the \code{graphlet_counts_1} argument is defined then
#' \code{graph_1} will not be used.
#'
#' @param graphlet_counts_compare Named list of pre-generated graphlet counts
#' for the remaining query graphs. If the \code{graphlet_counts_compare}
#' argument is defined then \code{graphs_compare} will not be used.
#'
#' @return Netdis statistics between graph_1 and graph_2 for graphlet sizes
#' up to and including max_graphlet_size
#' @export
netdis_one_to_many <- function(graph_1 = NULL,
                               graphs_compare = NULL,
                               ref_graph = 0,
                               max_graphlet_size = 4,
                               neighbourhood_size = 2,
                               min_ego_nodes = 3,
                               min_ego_edges = 1,
                               binning_fn = purrr::partial(
                                 binned_densities_adaptive,
                                 min_counts_per_interval = 5,
                                 num_intervals = 100
                               ),
                               bin_counts_fn = purrr::partial(
                                 density_binned_counts,
                                 agg_fn = mean,
                                 scale_fn = scale_graphlet_counts_ego
                               ),
                               exp_counts_fn = purrr::partial(
                                 netdis_expected_counts,
                                 scale_fn = count_graphlet_tuples
                               ),
                               graphlet_counts_1 = NULL,
                               graphlet_counts_compare = NULL) {
  ## ------------------------------------------------------------------------
  # Check arguments
  if (is.null(graph_1) & is.null(graphlet_counts_1)) {
    stop("One of graph_1 and graphlet_counts_1 must be supplied.")
  }
  if (is.null(graphs_compare) & is.null(graphlet_counts_compare)) {
    stop("One of graph_2 and graphlet_counts_2 must be supplied.")
  }
  
  ## ------------------------------------------------------------------------
  # Generate graphlet counts and bundle them into named list with format needed
  # for netdis_many_to_many.
  
  if (is.null(graphlet_counts_1)) {
    graphlet_counts_1 <- count_graphlets_ego(
      graph_1,
      max_graphlet_size = max_graphlet_size,
      neighbourhood_size = neighbourhood_size,
      min_ego_nodes = min_ego_nodes,
      min_ego_edges = min_ego_edges,
      return_ego_networks = FALSE
    )
  }
  rm(graph_1)
  
  if (is.null(graphlet_counts_compare)) {
    graphlet_counts_compare <- purrr::map(
      graphs_compare,
      count_graphlets_ego,
      max_graphlet_size = max_graphlet_size,
      neighbourhood_size = neighbourhood_size,
      min_ego_nodes = min_ego_nodes,
      min_ego_edges = min_ego_edges,
      return_ego_networks = FALSE
    )
  }
  rm(graphs_compare)
  
  graphlet_counts <- append(graphlet_counts_compare,
                            list(graph_1 = graphlet_counts_1),
                            after = 0
  )
  
  ## ------------------------------------------------------------------------
  # calculate netdis
  result <- netdis_many_to_many(
    graphs = NULL,
    ref_graph = ref_graph,
    comparisons = "one-to-many",
    max_graphlet_size = 4,
    neighbourhood_size = 2,
    min_ego_nodes = 3,
    min_ego_edges = 1,
    binning_fn = binning_fn,
    bin_counts_fn = bin_counts_fn,
    exp_counts_fn = exp_counts_fn,
    graphlet_counts = graphlet_counts
  )
  
  ## ------------------------------------------------------------------------
  # restructure netdis_many_to_many output
  colnames(result$netdis) <- result$comp_spec$name_b
  result$netdis
}


#' Netdis between all graph pairs
#'
#' @param graphs A named list of simplified igraph graph objects (undirected
#' graphs excluding loops, multiple edges and isolated vertices), such as those
#' obtained by using \code{read_simple_graphs}.
#'
#' @param ref_graph Controls how expected counts are calculated. Either:
#' 1) A numeric value - used as a constant expected counts value for all query
#' graphs (DEFAULT: 0).
#' 2) A simplified \code{igraph} object - used as a reference graph from which
#' expected counts are calculated for all query graphs.
#' 3) NULL - Expected counts will be calculated based on the properties of the
#' query graphs themselves.
#'
#' @param graphlet_counts Pre-generated graphlet counts. If the
#' \code{graphlet_counts} argument is defined then \code{graphs} will not be
#' used.
#' A named list of matrices containing counts of each graphlet (columns) for
#' each ego-network in the input graph (rows). Columns are labelled with
#' graphlet IDs and rows are labelled with the ID of the central node in each
#' ego-network. As well as graphlet counts, each matrix  must contain an
#' additional column labelled "N" including the node count for
#' each ego network.
#'
#' @param graphlet_counts_ref Pre-generated reference graphlet counts. If the
#' \code{graphlet_counts_ref} argument is defined then \code{ref_graph} will not
#' be used.
#' 
#' @param comparisons Which comparisons to perform between graphs.
#' Can be "many-to-many" (all pairwise combinations) or "one-to-many"
#' (compare first graph in graphs to all other graphs.)
#'
#' @param max_graphlet_size Generate graphlets up to this size.
#'
#' @param neighbourhood_size Ego network neighbourhood size.
#'
#' @param min_ego_nodes Filter ego networks which have fewer
#' than min_ego_nodes nodes.
#'
#' @param min_ego_edges Filter ego networks which have fewer
#' than min_ego_edges edges.
#'
#' @param binning_fn Function used to bin ego network densities. Takes densities
#' as its single argument, and returns a named list including keys \code{breaks}
#' (list of bin edges) and \code{interval_indexes} (density bin index for each
#' ego network). (Default: \code{binned_densities_adaptive} with
#' \code{min_counts_per_interval = 5} and \code{num_intervals = 100}).
#'
#' @param bin_counts_fn Function used to calculate expected graphlet counts in
#' each density bin. Takes \code{graphlet_counts}, \code{interval_indexes}
#' (bin indexes) and \code{max_graphlet_size} as arguments.
#' (Default: \code{density_binned_counts} with \code{agg_fn = mean} and
#' \code{scale_fn = scale_graphlet_counts_ego}, which mirrors the
#' approach used in the original netdis paper).
#'
#' @param exp_counts_fn Function used to map from binned reference counts to
#' expected counts for each graphlet in each ego network of the query graphs.
#' Takes \code{ego_networks}, \code{density_bin_breaks},
#' \code{binned_graphlet_counts}, and \code{max_graphlet_size} as arguments.
#' (Default: \code{netdis_expected_counts} with
#' \code{scale_fn = count_graphlet_tuples}, which mirrors the approach used in
#' the original netdis paper).
#'
#' @return Netdis statistics between query graphs for graphlet sizes
#' up to and including max_graphlet_size.
#'
#' @export
netdis_many_to_many <- function(graphs = NULL,
                                ref_graph = 0,
                                comparisons = "many-to-many",
                                max_graphlet_size = 4,
                                neighbourhood_size = 2,
                                min_ego_nodes = 3,
                                min_ego_edges = 1,
                                binning_fn = purrr::partial(
                                  binned_densities_adaptive,
                                  min_counts_per_interval = 5,
                                  num_intervals = 100
                                ),
                                bin_counts_fn = purrr::partial(
                                  density_binned_counts,
                                  agg_fn = mean,
                                  scale_fn = scale_graphlet_counts_ego
                                ),
                                exp_counts_fn = purrr::partial(
                                  netdis_expected_counts,
                                  scale_fn = count_graphlet_tuples
                                ),
                                graphlet_counts = NULL,
                                graphlet_counts_ref = NULL) {
  
  ## ------------------------------------------------------------------------
  # Check arguments
  if (is.null(graphs) & is.null(graphlet_counts)) {
    stop("One of graphs and graphlet_counts must be supplied.")
  }
  
  ## ------------------------------------------------------------------------
  # Generate ego networks and count graphlets for query graphs.
  # But if graphlet counts have already been provided we can skip this step.
  if (is.null(graphlet_counts)) {
    graphlet_counts <- purrr::map(
      graphs,
      count_graphlets_ego,
      max_graphlet_size = max_graphlet_size,
      neighbourhood_size = neighbourhood_size,
      min_ego_nodes = min_ego_nodes,
      min_ego_edges = min_ego_edges,
      return_ego_networks = FALSE
    )
  }
  rm(graphs)
  
  ## ------------------------------------------------------------------------
  # Centre counts
  # If there are no graphlet_counts_ref, and a number has been passed as ref_graph, treat it as a constant expected
  # counts value (e.g. if ref_graph = 0 then no centring of counts).
  if (is.numeric(ref_graph) && length(ref_graph) == 1 && is.null(graphlet_counts_ref)) {
    centred_graphlet_counts <- purrr::map(
      graphlet_counts,
      netdis_centred_graphlet_counts,
      ref_ego_density_bins = NULL,
      ref_binned_graphlet_counts = ref_graph,
      binning_fn = NULL,
      bin_counts_fn = NULL,
      exp_counts_fn = NULL,
      max_graphlet_size = max_graphlet_size
    )
    
    ## ------------------------------------------------------------------------
    # If there are no graphlet_counts_ref, and If a reference graph passed, use it to calculate expected counts for all
    # query graphs.
  } else if (!is.null(ref_graph) | !is.null(graphlet_counts_ref)) {
    
    # Generate ego networks and calculate graphlet counts
    # But if some ref graphlet counts provided can skip this step
    if (is.null(graphlet_counts_ref)) {
      graphlet_counts_ref <- count_graphlets_ego(
        ref_graph,
        max_graphlet_size = max_graphlet_size,
        neighbourhood_size = neighbourhood_size,
        min_ego_nodes = min_ego_nodes,
        min_ego_edges = min_ego_edges,
        return_ego_networks = FALSE
      )
    }
    rm(ref_graph)
    
    # Get ego-network densities
    densities_ref <- ego_network_density(graphlet_counts_ref)
    
    # bin ref ego-network densities
    binned_densities <- binning_fn(densities_ref)
    
    ref_ego_density_bins <- binned_densities$breaks
    
    # Average ref graphlet counts across density bins
    ref_binned_graphlet_counts <- bin_counts_fn(
      graphlet_counts_ref,
      binned_densities$interval_indexes,
      max_graphlet_size = max_graphlet_size
    )
    
    # Calculate centred counts using ref graph
    centred_graphlet_counts <- purrr::map(
      graphlet_counts,
      netdis_centred_graphlet_counts,
      ref_ego_density_bins = ref_ego_density_bins,
      ref_binned_graphlet_counts = ref_binned_graphlet_counts,
      binning_fn = binning_fn,
      bin_counts_fn = bin_counts_fn,
      exp_counts_fn = exp_counts_fn,
      max_graphlet_size = max_graphlet_size
    )
    
    ## ------------------------------------------------------------------------
    # If no reference passed, calculate expected counts using query networks
    # themselves.
  } else {
    centred_graphlet_counts <- purrr::map(
      graphlet_counts,
      netdis_centred_graphlet_counts,
      ref_ego_density_bins = NULL,
      ref_binned_graphlet_counts = NULL,
      binning_fn = binning_fn,
      bin_counts_fn = bin_counts_fn,
      exp_counts_fn = exp_counts_fn,
      max_graphlet_size = max_graphlet_size
    )
  }
  rm(graphlet_counts)
  
  ## ------------------------------------------------------------------------
  # Sum centred graphlet counts across all ego networks
  sum_graphlet_counts <- lapply(centred_graphlet_counts, colSums)
  
  rm(centred_graphlet_counts)
  
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
    SIMPLIFY = TRUE
  )
  
  list(netdis = results, comp_spec = comp_spec)
}

#' Netdis - for one graphlet size
#'
#' Calculate Netdis statistic between two graphs from their Centred Graphlet
#' Counts (generated using \code{netdis_centred_graphlet_counts}) for graphlets
#' of size \code{graphlet_size}.
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

#' Netdis - for all graphlet sizes up to max_graphlet_size
#'
#' Calculate Netdis statistic between two graphs from their Centred Graphlet
#' Counts (generated using \code{netdis_centred_graphlet_counts}) for all
#' graphlet sizes up to \code{max_graphlet_size}.
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


#' netdis_centred_graphlet_counts
#'
#' Calculate expected graphlet counts for each ego network in a query graph and
#' centre the actual counts by subtracting those calculated expected count
#' values.
#' @param graphlet_counts Ego network graphlet counts for a query graph
#'
#' @param ref_ego_density_bins Either a list of previously calculated ego
#' network density bin edges from a reference network, or \code{NULL}, in
#' which case density bins are generated using the query graph itself.
#'
#' @param ref_binned_graphlet_counts Either expected graphlet counts for each
#' ego network density bin from a reference network (a matrix with columns
#' labelled by graphlet ID and rows by density bin index), \code{NULL}, in
#' which case density binned counts are generated using the query graph itself,
#' or a constant numeric value to subtract from all graphlet counts.
#'
#' @param binning_fn Function used to bin ego network densities. Only needed if
#' \code{ref_ego_density_bins} and \code{ref_binned_graphlet_counts} are
#' \code{NULL}. Takes densities as its single argument, and returns a named list
#' including keys \code{breaks} (list of bin edges) and \code{interval_indexes}
#' (density bin index for each ego network).
#'
#' @param bin_counts_fn Function used to calculate expected graphlet counts in
#' each density bin. Only needed if  \code{ref_ego_density_bins} and
#' \code{ref_binned_graphlet_counts} are \code{NULL}. Takes
#' \code{graphlet_counts}, \code{interval_indexes} (bin indexes) and
#' \code{max_graphlet_size} as arguments.
#'
#' @param exp_counts_fn Function used to map from binned reference counts to
#' expected counts for each graphlet in each ego network of the query graphs.
#' Takes \code{ego_networks}, \code{density_bin_breaks},
#' \code{binned_graphlet_counts}, and \code{max_graphlet_size} as arguments.
#'
#' @param max_graphlet_size max graphlet size to calculate centred counts for.
#'
#' @return graphlet_counts minus exp_graphlet_counts for graphlets up to size
#' max_graphlet_size.
#' @export
netdis_centred_graphlet_counts <- function(
  graphlet_counts,
  ref_ego_density_bins,
  ref_binned_graphlet_counts,
  binning_fn,
  bin_counts_fn,
  exp_counts_fn,
  max_graphlet_size) {
  
  ## ------------------------------------------------------------------------
  # If a number has been passed as ref_binned_graphlet_counts, treat it as a
  # constant expected counts value (e.g. if ref_binned_graphlet_counts = 0
  # then no centring of counts).
  if (is.numeric(ref_binned_graphlet_counts) &&
      length(ref_binned_graphlet_counts) == 1) {
    exp_graphlet_counts <- netdis_const_expected_counts(
      graphlet_counts,
      const = ref_binned_graphlet_counts
    )
    
    ## ------------------------------------------------------------------------
    # If reference bins and counts passed, use them to calculate
    # expected counts
  } else if (!is.null(ref_ego_density_bins) &&
             !is.null(ref_binned_graphlet_counts)) {
    # Calculate expected graphlet counts (using ref
    # graph ego network density bins)
    exp_graphlet_counts <- exp_counts_fn(
      graphlet_counts,
      ref_ego_density_bins,
      ref_binned_graphlet_counts,
      max_graphlet_size = max_graphlet_size
    )
    
    ## ------------------------------------------------------------------------
    # If NULL passed as ref bins and counts, calculate expected counts using
    # query network itself.
  } else if (is.null(ref_ego_density_bins) &&
             is.null(ref_binned_graphlet_counts)) {
    # Get ego-network densities
    densities <- ego_network_density(graphlet_counts)
    
    # bin ref ego-network densities
    binned_densities <- binning_fn(densities)
    
    # extract bin breaks and indexes from binning results
    ego_density_bin_breaks <- binned_densities$breaks
    ego_density_bin_indexes <- binned_densities$interval_indexes
    
    # Calculate expected counts in each bin
    binned_graphlet_counts <- bin_counts_fn(
      graphlet_counts,
      ego_density_bin_indexes,
      max_graphlet_size = max_graphlet_size
    )
    
    # Calculate expected graphlet counts for each ego network
    exp_graphlet_counts <- exp_counts_fn(
      graphlet_counts,
      ego_density_bin_breaks,
      binned_graphlet_counts,
      max_graphlet_size = max_graphlet_size
    )
    
    ## ------------------------------------------------------------------------
    # Invalid combination of ref_ego_density_bins and ref_binned_graphlet_counts
  } else {
    stop("Invalid combination of ref_ego_density_bins and
         ref_binned_graphlet_counts. Options are:
         - Both NULL: calculate expected counts using query network.
         - List of bin edges and matrix of binned counts: Reference graph values
         for calculating expected counts.
         - Constant numeric ref_binned_graphlet_counts: Use as constant expected
         counts value.")
  }
  
  ## ------------------------------------------------------------------------
  # Subtract expected counts from actual graphlet counts
  netdis_subtract_exp_counts(
    graphlet_counts,
    exp_graphlet_counts,
    max_graphlet_size
  )
}


#' netdis_subtract_exp_counts
#'
#' Subtract expected graphlet counts from actual graphlet counts.
#'
#' @param graphlet_counts Matrix of graphlet counts (columns) for a
#' nummber of ego networks (rows).
#' @param exp_graphlet_counts Matrix of expected graphlet counts (columns) for a
#' nummber of ego networks (rows).
#' @param max_graphlet_size Do the subtraction for graphlets up to this size.
#' @export
netdis_subtract_exp_counts <- function(
  graphlet_counts,
  exp_graphlet_counts,
  max_graphlet_size) {
  
  # select columns for graphlets up to size max_graphlet_size
  id <- graphlet_key(max_graphlet_size)$id
  graphlet_counts <- graphlet_counts[, id]
  exp_graphlet_counts <- exp_graphlet_counts[, id]
  
  # Subtract expected counts from actual graphlet counts
  graphlet_counts - exp_graphlet_counts
}

#' netdis_expected_counts
#'
#' Calculates expected graphlet counts for each ego network based on its density
#' and pre-calculated reference density bins and graphlet counts for each bin.
#'
#' @param graphlet_counts Matrix of graphlet and node counts (columns) for a
#' nummber of ego networks (rows).
#' @param density_breaks Density values defining bin edges.
#' @param density_binned_reference_counts Reference network graphlet counts for
#' each density bin.
#' @param max_graphlet_size Determines the maximum size of graphlets to count.
#' Only graphlets containing up to \code{max_graphlet_size} nodes are counted.
#' @param scale_fn Optional function to scale calculated expected counts, taking
#' \code{graphlet_counts} and \code{max_graphlet_size} as arguments,
#' and returning a scale factor that the looked up
#' \code{density_binned_reference_counts} values will be multiplied by.
#'
#' @export
netdis_expected_counts <- function(
  graphlet_counts,
  density_breaks,
  density_binned_reference_counts,
  max_graphlet_size,
  scale_fn = NULL) {
  
  
  # Map over query graph ego-networks, using reference graph statistics to
  # calculate expected graphlet counts for each ego-network.
  expected_graphlet_counts <- t(apply(
    graphlet_counts, 1, netdis_expected_counts_ego,
    max_graphlet_size = max_graphlet_size,
    density_breaks = density_breaks,
    density_binned_reference_counts = density_binned_reference_counts,
    scale_fn = scale_fn
  ))
  
  expected_graphlet_counts
}

#' netdis_expected_counts_ego
#' INTERNAL FUNCTION - Do not call directly
#'
#' Calculates expected graphlet counts for one ego network based on its density
#' and pre-calculated reference density bins and graphlet counts for each bin.
#'
#' @param graphlet_counts Node and graphlet counts for an ego network.
#' @param max_graphlet_size Determines the maximum size of graphlets to count.
#' Only graphlets containing up to \code{max_graphlet_size} nodes are counted.
#' @param density_breaks Density values defining bin edges.
#' @param density_binned_reference_counts Reference network graphlet counts for
#' each density bin.
#' @param scale_fn Optional function to scale calculated expected counts, taking
#' \code{graphlet_counts} and \code{max_graphlet_size} as arguments, and
#' returning a scale factor that the looked up
#' \code{density_binned_reference_counts} values will be multiplied by.
#'
netdis_expected_counts_ego <- function(graphlet_counts,
                                       max_graphlet_size,
                                       density_breaks,
                                       density_binned_reference_counts,
                                       scale_fn = NULL) {
  
  # Look up average scaled graphlet counts for graphs of similar density
  # in the reference graph
  query_density <- density_from_counts(graphlet_counts)
  matched_density_index <- interval_index(query_density, density_breaks)
  
  matched_reference_counts <-
    density_binned_reference_counts[matched_density_index, ]
  
  if (!is.null(scale_fn)) {
    # Scale reference counts e.g. by multiplying the
    # reference count for each graphlet by the number
    # of possible sets of k nodes in the query graph,
    # where k is the number of nodes in the graphlet.
    matched_reference_counts <- matched_reference_counts *
      scale_fn(graphlet_counts, max_graphlet_size)
  }
  matched_reference_counts
}

#' mean_density_binned_graphlet_counts
#'
#' Calculate mean (dy default) graphlet counts for ego networks in each density
#' bin.
#'
#' @param graphlet_counts Graphlet counts for a number of ego_networks.
#' @param density_interval_indexes Density bin index for
#' each ego network in graphlet_counts.
#' @param agg_fn Function to aggregate counts in each bin
#' (default \code{agg_fn = mean}).
#'
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
#' the full range of possible density values (0 to 1).
#' @param densities Ego network density values (only used to return
#' a list of indexes of the required length.)
#' @export
single_density_bin <- function(densities) {
  list(
    densities = densities,
    interval_indexes = rep(1, length(densities)),
    breaks = c(0, 1)
  )
}

#' Used to calculate aggregated graphlet counts for each density bin.
#'
#' @param graphlet_counts Graphlet and node counts (columns) for a number of
#' ego_networks (rows).
#' @param density_interval_indexes Density bin index for
#' each ego network.
#' @param agg_fn Function to aggregate counts in each bin
#' (default \code{agg_fn = mean}).
#' @param scale_fn Optional function to apply a transformation/scaling
#' to the raw graphlet_counts. Must have arguments \code{graphlet_counts} and
#' \code{max_graphlet_size}, and return a transformed \code{graphlet_counts}
#' object with the same number of rows as the input, and columns for all
#' graphlets up to \code{max_graphlet_size}.
#' @param max_graphlet_size Optionally passed and used by scale_fn.
#'
#' @export
density_binned_counts <- function(graphlet_counts,
                                  density_interval_indexes,
                                  agg_fn = mean,
                                  scale_fn = NULL,
                                  max_graphlet_size = NULL) {
  if (!is.null(scale_fn)) {
    # Scale ego-network graphlet counts e.g.
    # by dividing by total number of k-tuples in
    # ego-network (where k is graphlet size)
    graphlet_counts <- scale_fn(graphlet_counts,
                                max_graphlet_size = max_graphlet_size
    )
  }
  
  mean_density_binned_graphlet_counts(graphlet_counts,
                                      density_interval_indexes,
                                      agg_fn = agg_fn
  )
}

#' INTERNAL FUNCTION - DO NOT CALL DIRECTLY
#' Used by \code{density_binned_counts_gp}
#' Calculate expected counts with geometric poisson (Polya-Aeppli)
#' approximation for a single density bin.
#' @param bin_idx Density bin index to calculate expected counts for.
#' @param graphlet_counts Graphlet counts for a number of ego_networks.
#' @param density_interval_indexes Density bin indexes for each ego network in
#' \code{graphlet_counts}.
#' @param max_graphlet_size Determines the maximum size of graphlets
#' included in graphlet_counts.
exp_counts_bin_gp <- function(bin_idx, graphlet_counts,
                              density_interval_indexes,
                              max_graphlet_size) {
  # extract ego networks belonging to input density bin index
  counts <- graphlet_counts[density_interval_indexes == bin_idx, ]
  
  # mean graphlet counts in this density bin
  means <- colMeans(counts)
  
  # subtract mean graphlet counts from actual graphlet counts
  mean_sub_counts <- sweep(counts, 2, means)
  
  # variance in graphlet counts across ego networks in this density bin
  Vd_sq <- colSums(mean_sub_counts^2) / (nrow(mean_sub_counts) - 1)
  
  # GP theta parameter for each graphlet id in this density bin
  theta_d <- 2 * means / (Vd_sq + means)
  
  exp_counts_dk <- vector()
  for (k in 2:max_graphlet_size) {
    graphlet_idx <- graphlet_ids_for_size(k)
    
    # GP lambda parameter for graphlet size k in this density bin
    lambda_dk <- mean(2 * means[graphlet_idx]^2 /
                        (Vd_sq[graphlet_idx] + means[graphlet_idx]),
                      na.rm = TRUE
    )
    
    # Expected counts for graphlet size k in this density bin
    exp_counts_dk <- append(
      exp_counts_dk,
      lambda_dk / theta_d[graphlet_idx]
    )
  }
  
  exp_counts_dk
}

#' Calculate expected counts in density bins using the
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
  
  # calculate expected counts for each density bin index
  nbins <- length(unique(density_interval_indexes))
  expected_counts_bin <- t(sapply(
    1:nbins,
    exp_counts_bin_gp,
    graphlet_counts = graphlet_counts,
    density_interval_indexes = density_interval_indexes,
    max_graphlet_size = max_graphlet_size
  ))
  
  # remove NAs caused by bins with zero counts for a graphlet
  expected_counts_bin[is.nan(expected_counts_bin)] <- 0
  
  expected_counts_bin
}

#' Create matrix of constant value to use as expected counts.
#' @param graphlet_counts Ego network graphlet counts matrix to create expected
#' counts for.
#' @param const Constant expected counts value to use.
#' @return Counts of value const with same shape and names as graphlet_counts.
netdis_const_expected_counts <- function(graphlet_counts, const) {
  exp_counts <- graphlet_counts
  exp_counts[, ] <- const
  exp_counts
}


#' Replace zero values in a vector with ones. Used by
#' \code{scale_graphlet_count} to prevent divide by
#' zero errors.
#' @param v A vector.
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
  graphlet_count[, colnames(graphlet_tuples)] / zeros_to_ones(graphlet_tuples)
}


#' Run count_graphlet_tuples across pre-computed ego networks.
#' @param graphlet_counts Matrix of graphlet and node counts (columns) for a
#' number of ego networks (rows).
#' @param max_graphlet_size Determines the maximum size of graphlets included
#' in the tuple counts.
#' @export
count_graphlet_tuples_ego <- function(graphlet_counts, max_graphlet_size) {
  graphlet_tuple_counts <-
    t(apply(graphlet_counts, 1,
            count_graphlet_tuples,
            max_graphlet_size = max_graphlet_size
    ))
  
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
#' @param graph_graphlet_counts Node and graphlet counts for a single graph.
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
}
