

graph_To_dhist <- function(G,max_graphlet_size=4,feature_type=')orbit')
{
    counts1 <- graph_to_graphletCount(G,feature_type,max_graphlet_size,ego_neighbourhood_size)
    dhist1 <- graph_features_to_histograms(counts1)

}

graph_Pair_To_NetEMD <- function(G1,G2,feature_type = 'orbit',max_graphlet_size = 4)
{
    ## Construct histograms for each of the counts
    counts1 <- graph_to_graphletCount(G1,feature_type,max_graphlet_size,ego_neighbourhood_size)
    counts2 <- graph_to_graphletCount(G2,feature_type,max_graphlet_size,ego_neighbourhood_size)
    ## Make discrete histograms (dhists) from the counts 
    dhist1 <- graph_features_to_histograms(counts1)
    dhist2 <- graph_features_to_histograms(counts2)
    ## Compute the netemd on these dhists
    dhist_Pair_To_NetEMD(dhist1,dhist2)
}

multiple_graphs_To_NetEMD <- function(graphs,feature_type = 'orbit',max_graphlet_size = 4)
{
    dhists <- multiple_graphs_To_dhists(graphs)
    ## Compute the netemd on these dhists
    multiple_dhists_To_netEMD(dhist1,dhist2)
}

features_Pair_To_NetEMD <- function(features1,features2,feature_type = 'orbit',max_graphlet_size = 4)
{
    ## Make discrete histograms (dhists) from the features
    dhist1 <- graph_features_to_histograms(features1)
    dhist2 <- graph_features_to_histograms(features2)
    ## Compute the netemd on these dhists
    dhist_Pair_To_NetEMD(dhist1,dhist2)
}

#' Graphlet-based counts
#' 
#' Generates graphlet-based degree distributions from \code{igraph} graph object,
#' using the ORCA fast graphlet orbit counting package.
#' @param graph A connected, undirected, simple graph as an \code{igraph} object. 
#' @param feature_type Type of graphlet-based feature to count: "graphlet"
#' counts the number of graphlets each node participates in; "orbit" calculates
#' the number of graphlet orbits each node participates in.
#' @param max_graphlet_size Determines the maximum size of graphlets to count. 
#' Only graphlets containing up to \code{max_graphlet_size} nodes will be counted.
#' @param ego_neighbourhood_size The number of steps from the source node to include
#' nodes for each ego-network.
#' @return List of graphlet-based degree distributions, with each distribution
#' represented as a \code{dhist} discrete histogram object.
#' @export
graph_to_graphletCount<- function(graph, feature_type = 'orbit', max_graphlet_size = 4, 
                ego_neighbourhood_size = 0){
  graph <- simplify_graph(graph)
  if(ego_neighbourhood_size > 0) {
    if(feature_type != 'graphlet') {
      stop("Feature type not supported for ego-networks")
    } else {
      out <- count_graphlets_ego(graph, max_graphlet_size = max_graphlet_size, 
                                 neighbourhood_size = ego_neighbourhood_size)
    }
  } else  if(feature_type == "orbit") {
    out <- count_orbits_per_node(graph, max_graphlet_size = max_graphlet_size)
  } else if(feature_type == "graphlet") {
    out <- count_graphlets_per_node(graph, max_graphlet_size = max_graphlet_size)
  }
  else {
    stop('gdd: unrecognised feature_type')
  }
  graph_features_to_histograms(out)
}

#' NetEMDs between all graph pairs using provided Graphlet-based Degree 
#' Distributions 
#' @param gdds List containing sets of Graphlet-based Degree Distributions for 
#' all graphs being compared
#' @param method The method to use to find the minimum EMD across all potential 
#' offsets for each pair of histograms. Default is "optimise" to use
#' R's built-in \code{stats::optimise} method to efficiently find the offset 
#' with the minimal EMD. However, this is not guaranteed to find the global 
#' minimum if multiple local minima EMDs exist. You can alternatively specify the 
#' "exhaustive" method, which will exhaustively evaluate the EMD between the 
#' histograms at all offsets that are candidates for the minimal EMD.
#' @param return_details Logical indicating whether to return the individual
#' minimal EMDs and associated offsets for all pairs of histograms
#' @param smoothing_window_width Width of "top-hat" smoothing window to apply to
#' "smear" point masses across a finite width in the real domain. Default is 0, 
#' which  results in no smoothing. Care should be taken to select a 
#' \code{smoothing_window_width} that is appropriate for the discrete domain 
#' (e.g.for the integer domain a width of 1 is the natural choice)
#' @param  mc.cores Number of cores to use for parallel processing. Defaults to
#' the \code{mc.cores} option set in the R environment.
#' @return NetEMD measures between all pairs of graphs for which GDDs 
#' were provided. Format of returned data depends on the \code{return_details}
#' parameter. If set to FALSE, a list is returned with the following named
#' elements:\code{net_emd}: a vector of NetEMDs for each pair of graphs, 
#' \code{comp_spec}: a comaprison specification table containing the graph names 
#' and indices within the input GDD list for each pair of graphs compared.
#' If \code{return_details} is set to FALSE, the list also contains the following
#' matrices for each graph pair: \code{min_emds}: the minimal EMD for each GDD 
#' used to compute the NetEMD, \code{min_offsets}: the associated offsets giving
#' the minimal EMD for each GDD
#' @export
multiple_dhists_To_netEMD <- function(
  gdds, method = "optimise", smoothing_window_width = 0, 
  return_details = FALSE, mc.cores = getOption("mc.cores", 2L)) {
  comp_spec <- cross_comparison_spec(gdds)
  # NOTE: mcapply only works on unix-like systems with system level forking 
  # capability. This means it will work on Linux and OSX, but not Windows.
  # For now, we just revert to single threaded operation on Windows
  # TODO: Look into using the parLappy function on Windows
  if(.Platform$OS.type != "unix") {
    # Force cores to 1 if system is not unix-like as it will not support 
    # forking
    mc.cores = 1
  }
  num_features <- length(gdds[[1]])
  out <- purrr::simplify(parallel::mcmapply(function(index_a, index_b) {net_emd(
    gdds[[index_a]], gdds[[index_b]], method = method, return_details = return_details,
    smoothing_window_width = smoothing_window_width)
    }, comp_spec$index_a, comp_spec$index_b, SIMPLIFY = FALSE, mc.cores = mc.cores))
  if(return_details) {
    net_emds <- purrr::simplify(purrr::map(out, ~.$net_emd))
    min_emds <- matrix(purrr::simplify(purrr::map(out, ~.$min_emds)), ncol = num_features, byrow = TRUE)
    colnames(min_emds) <- purrr::simplify(purrr::map(1:num_features, ~paste("MinEMD_O", .-1, sep = "")))
    min_offsets <- matrix(purrr::simplify(purrr::map(out, ~.$min_offsets)), ncol = num_features, byrow = TRUE)
    colnames(min_offsets) <- purrr::simplify(purrr::map(1:num_features, ~paste("MinOffsets_O", .-1, sep = "")))
    ret <- list(net_emds = net_emds, comp_spec = comp_spec, min_emds = min_emds, min_offsets = min_offsets)
  } else {
    net_emds <- out
    ret <- list(net_emds = net_emds, comp_spec = comp_spec)
  }
}


#' NetEMD Network Earth Mover's Distance
#' 
#' Calculates the mean minimum Earth Mover's Distance (EMD) between two sets of
#' discrete histograms after normalising each histogram to unit mass and variance.
#' This is calculated as follows:
#'   1. Normalise each histogram to have unit mass and unit variance
#'   2. Find the minimum EMD between each pair of histograms
#'   3. Take the average minimum EMD across all histogram pairs
#' @param dhists1 A \code{dhist} discrete histogram object or a list of such objects
#' @param dhists2 A \code{dhist} discrete histogram object or a list of such objects
#' @param method The method to use to find the minimum EMD across all potential 
#' offsets for each pair of histograms. Default is "optimise" to use
#' R's built-in \code{stats::optimise} method to efficiently find the offset 
#' with the minimal EMD. However, this is not guaranteed to find the global 
#' minimum if multiple local minima EMDs exist. You can alternatively specify the 
#' "exhaustive" method, which will exhaustively evaluate the EMD between the 
#' histograms at all offsets that are candidates for the minimal EMD.
#' @param return_details Logical indicating whether to return the individual
#' minimal EMDs and associated offsets for all pairs of histograms
#' @param smoothing_window_width Width of "top-hat" smoothing window to apply to
#' "smear" point masses across a finite width in the real domain. Default is 0, 
#' which  results in no smoothing. Care should be taken to select a 
#' \code{smoothing_window_width} that is appropriate for the discrete domain 
#' (e.g.for the integer domain a width of 1 is the natural choice)
#' @return NetEMD measure for the two sets of discrete histograms 
#' (\code{return_details = FALSE}) or a list with the following named elements
#' \code{net_emd}: the NetEMD for the set of histogram pairs, \code{min_emds}:  
#' the minimal EMD for each pair of histograms, \code{min_offsets}: the associated
#' offsets giving the minimal EMD for each pair of histograms
#' @export
dhist_Pair_To_NetEMD <- function(dhists1, dhists2, method = "optimise", 
                    return_details = FALSE, smoothing_window_width = 0) {
  # Require either a pair of "dhist" discrete histograms or two lists of "dhist"
  # discrete histograms
  pair_of_dhist_lists <- all(purrr::map_lgl(dhists1, is_dhist)) && all(purrr::map_lgl(dhists2, is_dhist))
  
  # If input is two lists of "dhist" discrete histograms, determine the minimum
  # EMD and associated offset for pairs of histograms taken from the same 
  # position in each list
  if(pair_of_dhist_lists) {
    details <- purrr::map2(dhists1, dhists2, function(dhist1, dhist2) {
      net_emd_single_pair(dhist1, dhist2, method = method,
                          smoothing_window_width = smoothing_window_width)
      })
    # Collect the minimum EMDs and associated offsets for all histogram pairs
    min_emds <- purrr::simplify(purrr::transpose(details)$min_emd)
    min_offsets <- purrr::simplify(purrr::transpose(details)$min_offset)
    # The NetEMD is the arithmetic mean of the minimum EMDs for each pair of 
    # histograms
    arithmetic_mean <- sum(min_emds) / length(min_emds)
    net_emd <- arithmetic_mean
    # Return just the NetEMD or a list including the NetEMD plus the details of
    # the minumum EMD and associated offsets for the individual histograms
    # Note that the offsets represent shifts after the histograms have been
    # scaled to unit variance
    if(return_details) {
      return(list(net_emd = net_emd, min_emds = min_emds, min_offsets = min_offsets))
    } else {
      return(arithmetic_mean)
    }
  }
  else {
    # Wrap each member of a single pair of histograms is a list and recursively
    # call this net_emd function. This ensures they are treated the same.
    return(dhist_Pair_To_NetEMD (list(dhists1), list(dhists2), method = method, 
                   return_details = return_details,
                   smoothing_window_width = smoothing_window_width))
  }
}

#' Load all graphs in a directory and calculates their Graphlet-based Degree
#' Distributions (GDDs)
#' 
#' Loads graphs from all files matching the given pattern in the given directory,
#' converts them to indexed edge lists compatible with the ORCA fast orbit 
#' counting package and calculates the specified set of graphlet-based degree 
#' distributions usingthe ORCA package.
#' @param source_dir Path to graph directory
#' @param format Format of graph files
#' @param pattern Filename pattern to match graph files
#' @param feature_type Type of graphlet-based degree distributions. Can be 
#' \code{graphlet} to count graphlets or \code{orbit} to count orbits.
#' @return A named list where each element contains a set of GDDs for a single
#' @param max_graphlet_size Maximum size of graphlets to use when generating GDD
#' @param ego_neighbourhood_size The number of steps from the source node to 
#' include nodes for each ego-network. If set to 0, ego-networks will not be 
#' used
#' @param  mc.cores Number of cores to use for parallel processing. Defaults to
#' the \code{mc.cores} option set in the R environment.
#' @return A named list where each element contains a set of GDDs for a single
#' graph from the source directory. Each set of GDDs is itself a named list,  
#' where each GDD element is a \code{dhist} discrete histogram object.
#' @export
graphFolder_To_dhists <- function(
  source_dir, format = "ncol", pattern = ".txt", feature_type = "orbit", 
  max_graphlet_size = 4, ego_neighbourhood_size = 0,
  mc.cores = getOption("mc.cores", 2L)) {
  # Create function to read graph from file and generate GDD
  graphs <- read_simple_graphs( source_dir = source_dir, format = format,
                               pattern = pattern)
    multiple_graphs_To_dhists(graphs, feature_type = feature_type,
                            max_graphlet_size = max_graphlet_size,
                            ego_neighbourhood_size = ego_neighbourhood_size,
                            mc.cores = mc.cores)
}

multiple_graphs_To_dhists <- function(graphs,
  feature_type = "orbit", 
  max_graphlet_size = 4, ego_neighbourhood_size = 0,
  mc.cores = getOption("mc.cores", 2L)) {
  # Calculate specified GDDs for each graph
  # NOTE: mcapply only works on unix-like systems with system level forking 
  # capability. This means it will work on Linux and OSX, but not Windows.
  # For now, we just revert to single threaded operation on Windows
  # TODO: Look into using the parLappy function on Windows
  if(.Platform$OS.type != "unix") {
    # Force cores to 1 if system is not unix-like as it will not support 
    # forking
    mc.cores = 1
  }
  parallel::mcmapply(graph_To_dhist, graphs, MoreArgs = 
                       list(feature_type = feature_type, 
                            max_graphlet_size = max_graphlet_size,
                            ego_neighbourhood_size = ego_neighbourhood_size), 
                     SIMPLIFY = FALSE, mc.cores = mc.cores)
}


multiple_features_To_dhists <- function(features,
  mc.cores = getOption("mc.cores", 2L)) {
  # Calculate specified GDDs for each graph
  # NOTE: mcapply only works on unix-like systems with system level forking 
  # capability. This means it will work on Linux and OSX, but not Windows.
  # For now, we just revert to single threaded operation on Windows
  # TODO: Look into using the parLappy function on Windows
  if(.Platform$OS.type != "unix") {
    # Force cores to 1 if system is not unix-like as it will not support 
    # forking
    mc.cores = 1
  }
  parallel::mcmapply(graph_features_to_histograms, features,mc.cores = mc.cores)
}
