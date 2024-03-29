#' NetEMD Network Earth Mover's Distance between a pair of networks.
#'
#' Calculates the network Earth Mover's Distance (EMD) between 
#' two sets of network features. This is done by individually normalising the distribution
#' of each feature so that they have unit mass and unit variance. Then the minimun EMD between the same pair of features (one for each corresponding graph) is calculated by considering all possible translations of the feature distributions. Finally the average over all features is reported.
#' This is calculated as follows:
#'   1. Normalise each feature histogram to have unit mass and unit variance.
#'   2. For each feature, find the minimum EMD between each pair of histograms considering all possible  histogram translations.
#'   3. Take the average minimum EMD across all features.
#' @param graph_1 A network/graph object from the \code{igraph} package. \code{graph_1} can be set to \code{NULL} (default) if \code{dhists_1} is provided.
#' @param graph_2 A network/graph object from the \code{igraph} package. \code{graph_2} can be set to \code{NULL} (default) if \code{dhists_2} is provided.
#' @param dhists_1 Either, a \code{dhist} discrete histogram object, or list of such objects, or a matrix of network features (each column representing a feature). \code{dhists_1} can be set to \code{NULL} (default) if \code{graph_1} is provided. A \code{dhist} object can be obtained from \code{graph_features_to_histograms}.
#' @param dhists_2 Same as \code{dhists_1}.
#' @param method The method to be used to find the minimum EMD across all potential
#' offsets for each pair of histograms. Default is "optimise" to use
#' R's built-in \code{stats::optimise} method to efficiently find the offset
#' with the minimal EMD. However, this is not guaranteed to find the global
#' minimum if multiple local minima EMDs exist. You can alternatively specify the
#' "exhaustive" method, which will exhaustively evaluate the EMD between the
#' histograms at all offsets that are candidates for the minimal EMD at the cost of computational time.
#' @param return_details Logical indicating whether to return the individual
#' minimal EMDs and associated offsets for all pairs of histograms.
#' @param smoothing_window_width Width of "top-hat" smoothing window to apply to
#' "smear" point masses across a finite width in the real domain. Default is 0,
#' which  results in no smoothing. Care should be taken to select a
#' \code{smoothing_window_width} that is appropriate for the discrete domain
#' (e.g.for the integer domain a width of 1 is the natural choice).
#' @param feature_type Type of graphlet-based feature to count: "graphlet"
#' counts the number of graphlets each node participates in; "orbit" (default) calculates
#' the number of graphlet orbits each node participates in.
#' @param max_graphlet_size Determines the maximum size of graphlets to count.
#' Only graphlets containing up to \code{max_graphlet_size} nodes will be
#' counted. Possible values are 4, and 5 (default).
#' @param ego_neighbourhood_size The number of steps from the source node to
#' include nodes for each ego-network. NetEmd was proposed for individual nodes alone, hence the default value is 0.
#' @return NetEMD measure for the two sets of discrete histograms (or graphs). If
#' (\code{return_details = FALSE}) then a list with the following named elements is returned
#' \code{net_emd}: the NetEMD for the set of histogram pairs (or graphs), \code{min_emds}:
#' the minimal EMD for each pair of histograms, \code{min_offsets}: the associated
#' offsets giving the minimal EMD for each pair of histograms
#' @examples
#'  require(igraph)
#'  graph_1 <- graph.lattice(c(8,8)) 
#'  graph_2 <- graph.lattice(c(44,44)) 
#'  netemd_one_to_one(graph_1=graph_1,graph_2=graph_2,feature_type="orbit",max_graphlet_size=5)
#'  
#'  #Providing a matrix of network features
#'  props_a= count_orbits_per_node(graph = graph_1,max_graphlet_size = 5)
#'  props_b= count_orbits_per_node(graph = graph_2,max_graphlet_size = 5)
#'  
#'  netemd_one_to_one(dhists_1=props_a, dhists_2=props_b,smoothing_window_width = 1)
#'  
#'  #Providing the network features as lists of dhist objects
#'  dhists_1<- graph_features_to_histograms(props_a)
#'  dhists_2<- graph_features_to_histograms(props_b)
#'  
#'  netemd_one_to_one(dhists_1=dhists_1, dhists_2=dhists_2)
#'  
#'  
#'  # A variation of NetEmd: Using the Laplacian spectrum 
#'  #Laplacian
#'  Lapg_1 <- igraph::laplacian_matrix(graph = graph_1,normalized = FALSE,sparse = FALSE)
#'  Lapg_2 <- igraph::laplacian_matrix(graph = graph_2,normalized = FALSE,sparse = FALSE)
#'  
#'  #Normalized Laplacian
#'  NLapg_1 <- igraph::laplacian_matrix(graph = graph_1,normalized = TRUE,sparse = FALSE)
#'  NLapg_2 <- igraph::laplacian_matrix(graph = graph_2,normalized = TRUE,sparse = FALSE)
#'  
#'  #Spectra (This may take a couple of minutes).
#'  props_1 <- cbind(L.Spectra= eigen(Lapg_1)$values, NL.Spectra= eigen(NLapg_1)$values) 
#'  props_2 <- cbind(L.Spectra= eigen(Lapg_2)$values, NL.Spectra= eigen(NLapg_2)$values) 
#'  
#'  netemd_one_to_one(dhists_1 = props_1,dhists_2 = props_2,smoothing_window_width = 0)#Use of smoothing window 1 is given for discrete integer distributions. If the network features are considered continuous variables smoothing_window_width equal to zero is recommended.
#'  
#' @export
netemd_one_to_one <- function(graph_1=NULL,graph_2=NULL,dhists_1=NULL, dhists_2=NULL, method = "optimise",
                              return_details = FALSE, smoothing_window_width = 0,feature_type="orbit",max_graphlet_size = 5,ego_neighbourhood_size = 0) {
  ## ------------------------------------------------------------------------
  # Check arguments 1
  if (!igraph::is.igraph(graph_1) & is.null(dhists_1)) {
    stop("One of graph_1 or dhists_1 must be supplied.")
  }
  if (!igraph::is.igraph(graph_2) & is.null(dhists_2)) {
    stop("One of graph_2 or dhists_2 must be supplied.")
  }
  ## ------------------------------------------------------------------------
  # Check arguments 2
  # If dhists_1 is a matrix of network features then transform them to dhist objects.
  if(is.matrix(dhists_1)){
    dhists_1 <- graph_features_to_histograms(dhists_1)
  }
  if(is.matrix(dhists_2)){
    dhists_2 <- graph_features_to_histograms(dhists_2)
  }
  ## ------------------------------------------------------------------------
  # Check arguments 3
  #If input is graph then get graphlet counts
  if(igraph::is.igraph(graph_1)){
    if(!is.null(dhists_1)){warning("dhists_1 will be calculated from graph_1.")}
    dhists_1 <- gdd(graph = graph_1, feature_type = feature_type,
                    max_graphlet_size = max_graphlet_size,
                    ego_neighbourhood_size = ego_neighbourhood_size
    )
  }
  if(igraph::is.igraph(graph_2)){
    if(!is.null(dhists_2)){warning("dhists_2 will be calculated from graph_2.")}
    dhists_2 <- gdd(graph = graph_2, feature_type = feature_type,
                    max_graphlet_size = max_graphlet_size,
                    ego_neighbourhood_size = ego_neighbourhood_size
    )
  }
  
  rm(graph_1,graph_2)
  ## ------------------------------------------------------------------------
  # Require either a pair of "dhist" discrete histograms or two lists of "dhist"
  # discrete histograms
  pair_of_dhist_lists <- all(purrr::map_lgl(dhists_1, is_dhist)) && all(purrr::map_lgl(dhists_2, is_dhist))
  
  # If input is two lists of "dhist" discrete histograms, determine the minimum
  # EMD and associated offset for pairs of histograms taken from the same
  # position in each list
  if (pair_of_dhist_lists) {
    details <- purrr::map2(dhists_1, dhists_2, function(dhist1, dhist2) {
      netemd_single_pair(dhist1, dhist2,
                         method = method,
                         smoothing_window_width = smoothing_window_width
      )
    })
    # Collect the minimum EMDs and associated offsets for all histogram pairs
    min_emds <- purrr::simplify(purrr::transpose(details)$min_emd)
    min_offsets <- purrr::simplify(purrr::transpose(details)$min_offset)
    min_offsets_std <- purrr::simplify(purrr::transpose(details)$min_offset_std)
    # The NetEMD is the arithmetic mean of the minimum EMDs for each pair of
    # histograms
    arithmetic_mean <- sum(min_emds) / length(min_emds)
    net_emd <- arithmetic_mean
    # Return just the NetEMD or a list including the NetEMD plus the details of
    # the minumum EMD and associated offsets for the individual histograms
    # Note that the offsets represent shifts after the histograms have been
    # scaled to unit variance
    if (return_details) {
      return(list(net_emd = net_emd, min_emds = min_emds, min_offsets = min_offsets, min_offsets_std = min_offsets_std))
    } else {
      return(arithmetic_mean)
    }
  }
  else {
    # Wrap each member of a single pair of histograms is a list and recursively
    # call this net_emd function. This ensures they are treated the same.
    return(netemd_one_to_one(dhists_1 = list(dhists_1), dhists_2 = list(dhists_2),
                             method = method,
                             return_details = return_details,
                             smoothing_window_width = smoothing_window_width
    ))
  }
}


#' NetEMDs between all graph pairs using provided Graphlet-based Degree
#' Distributions
#' @param graphs A list of network/graph objects from the \code{igraph} package. \code{graphs} can be set to \code{NULL} (default) if \code{dhists} is provided.
#' @param dhists A list whose elements contain either: A  list of \code{dhist} discrete histogram objects for each graph, or a list  a matrix of network features (each column representing a feature). \code{dhists} can be set to \code{NULL} (default) if \code{graphs} is provided.  A \code{dhist} object can be obtained from \code{graph_features_to_histograms}.
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
#' (e.g.for the integer domain a width of 1 is the natural choice).
#' @param  mc.cores Number of cores to use for parallel processing. Defaults to
#' the \code{mc.cores} option set in the R environment.
#' @param feature_type Type of graphlet-based feature to count: "graphlet"
#' counts the number of graphlets each node participates in; "orbit" (default) calculates
#' the number of graphlet orbits each node participates in.
#' @param max_graphlet_size Determines the maximum size of graphlets to count.
#' Only graphlets containing up to \code{max_graphlet_size} nodes will be
#' counted. Possible values are 4, and 5 (default).
#' @param ego_neighbourhood_size The number of steps from the source node to
#' include nodes for each ego-network. NetEmd was proposed for individual nodes alone, hence the default value is 0.
#' @return NetEMD measures between all pairs of graphs for which features
#' were provided. Format of returned data depends on the \code{return_details}
#' parameter. If set to FALSE, a list is returned with the following named
#' elements:\code{net_emd}: a vector of NetEMDs for each pair of graphs,
#' \code{comp_spec}: a comparison specification table containing the graph names
#' and indices within the input GDD list for each pair of graphs compared.
#' If \code{return_details} is set to FALSE, the list also contains the following
#' matrices for each graph pair: \code{min_emds}: the minimal EMD for each GDD
#' used to compute the NetEMD, \code{min_offsets}: the associated offsets giving
#' the minimal EMD for each GDD
#' @export
netemd_many_to_many<- function(graphs=NULL,dhists=NULL, method = "optimise", smoothing_window_width = 0,
                               return_details = FALSE, mc.cores = getOption("mc.cores", 2L),feature_type="orbit",max_graphlet_size = 5,ego_neighbourhood_size = 0) {
  if(max_graphlet_size > 4 & mc.cores > 1) print(paste("This function will compute orbits of graphlets up to size 5 using ", mc.cores," cores. Depending on the density and size of the graphs, this may lead to a large compsumption of RAM."))
  
  # NOTE: mcapply only works on unix-like systems with system level forking
  # capability. This means it will work on Linux and OSX, but not Windows.
  # For now, we just revert to single threaded operation on Windows
  # TODO: Look into using the parLappy function on Windows
  if (.Platform$OS.type != "unix") {
    # Force cores to 1 if system is not unix-like as it will not support
    # forking
    mc.cores <- 1
  }
  ## ------------------------------------------------------------------------
  # Check arguments 1
  which_imput_type <- NULL
  if(!is.null(graphs) & is.null(dhists)){
    if ( !all(( unlist(sapply(X = graphs, FUN = igraph::is.igraph)) ) )  ) {
      stop("Graphs need to be igraph graph objects, or a list of dhists network features should be supplied.")
    }
    which_imput_type <- "Graphs"
  }
  if (!is.null(dhists) ) {
    if (all(( unlist(sapply(X = dhists, FUN = is.matrix)) ) )  ) {
      which_imput_type <- "Matrix"
    } 
    if ( all(( unlist(sapply(X = dhists, FUN = 
                             function(l){ all(( unlist(sapply(X = l, FUN = is_dhist)) ) ) }
    )) ) )  ) {
      which_imput_type <- "dhist"
    }
    if(is.null(which_imput_type)){
      warning("dhists does not conform to a Matrix or dhist class for all elmenents/netwroks in the list.")
    }
  }
  ## ------------------------------------------------------------------------
  # Check arguments 2
  # If dhists is a list of matrices of network features then transform them to dhist objects.
  if(which_imput_type == "Matrix"){
    dhists <- sapply(X = dhists,FUN = graph_features_to_histograms, simplify = FALSE )
  }
  ## ------------------------------------------------------------------------
  # Check arguments 3
  #If input is graph then get graphlet counts
  if(which_imput_type == "Graphs"){
    dhists <- parallel::mcmapply(gdd, graphs,
                                 MoreArgs =
                                   list(
                                     feature_type = feature_type,
                                     max_graphlet_size = max_graphlet_size,
                                     ego_neighbourhood_size = ego_neighbourhood_size
                                   ),
                                 SIMPLIFY = FALSE, mc.cores = mc.cores
    )
  }
  rm(graphs)
  ## ------------------------------------------------------------------------
  # Check arguments 4
  #cross_comparison_spec was coded to require names!
  if(is.null(names(dhists))){
    names(dhists) <- paste("Net",1:length(dhists),sep = "")
  }
  ## ------------------------------------------------------------------------
  comp_spec <- cross_comparison_spec(dhists)
  num_features <- length(dhists[[1]])
  out <- purrr::simplify(parallel::mcmapply(function(index_a, index_b) {
    netemd_one_to_one(dhists_1 =  dhists[[index_a]], dhists_2 =  dhists[[index_b]],
                      method = method, return_details = return_details,
                      smoothing_window_width = smoothing_window_width
    )
  }, comp_spec$index_a, comp_spec$index_b, SIMPLIFY = FALSE, mc.cores = mc.cores))
  if (return_details) {
    net_emds <- purrr::simplify(purrr::map(out, ~ .$net_emd))
    min_emds <- matrix(purrr::simplify(purrr::map(out, ~ .$min_emds)), ncol = num_features, byrow = TRUE)
    colnames(min_emds) <- purrr::simplify(purrr::map(1:num_features, ~ paste("MinEMD_O", . - 1, sep = "")))
    min_offsets <- matrix(purrr::simplify(purrr::map(out, ~ .$min_offsets)), ncol = num_features, byrow = TRUE)
    colnames(min_offsets) <- purrr::simplify(purrr::map(1:num_features, ~ paste("MinOffsets_O", . - 1, sep = "")))
    min_offsets_std <- matrix(purrr::simplify(purrr::map(out, ~ .$min_offsets_std)), ncol = num_features, byrow = TRUE)
    colnames(min_offsets_std) <- purrr::simplify(purrr::map(1:num_features, ~ paste("MinOffsetsStd_O", . - 1, sep = "")))
    ret <- list(netemds = net_emds, comp_spec = comp_spec, min_emds = min_emds, min_offsets = min_offsets, min_offsets_std = min_offsets_std)
  } else {
    net_emds <- out
    ret <- list(netemds = net_emds, comp_spec = comp_spec)
  }
  return(ret)
}

#' Internal function to compute the minimum Earth Mover's Distance between standarized and translated histograms
#'
#' Calculates the minimum Earth Mover's Distance (EMD) between two
#' discrete histograms after normalising each histogram to unit mass and variance.
#' This is calculated as follows:
#'   1. Normalise each histogram to have unit mass and unit variance
#'   2. Find the minimum EMD between the histograms
#' @param dhists_1 A \code{dhist} discrete histogram object or a list of such objects
#' @param dhists_2 A \code{dhist} discrete histogram object or a list of such objects
#' @param method The method to use to find the minimum EMD across all potential
#' offsets for each pair of histograms. Default is "optimise" to use
#' R's built-in \code{stats::optimise} method to efficiently find the offset
#' with the minimal EMD. However, this is not guaranteed to find the global
#' minimum if multiple local minima EMDs exist. You can alternatively specify the
#' "exhaustive" method, which will exhaustively evaluate the EMD between the
#' histograms at all offsets that are candidates for the minimal EMD.
#' @param smoothing_window_width Width of "top-hat" smoothing window to apply to
#' "smear" point masses across a finite width in the real domain. Default is 0,
#' which  results in no smoothing. Care should be taken to select a
#' \code{smoothing_window_width} that is appropriate for the discrete domain
#' (e.g.for the integer domain a width of 1 is the natural choice)
#' @return A list with the following named elements
#' \code{net_emd}: the NetEMD for the set of histogram pairs, \code{min_offsets}: the associated
#' offsets giving the minimal EMD for each pair of histograms and \code{min_offset_std}: Offset used in the standardised histograms.
#' @examples 
#'  require(igraph)
#'  goldstd_1 <- graph.lattice(c(8,8)) 
#'  goldstd_2 <- graph.lattice(c(44,44)) 
#'  props_1 <- count_orbits_per_node(graph = goldstd_1,max_graphlet_size = 5)
#'  props_2 <- count_orbits_per_node(graph = goldstd_2,max_graphlet_size = 5)
#'  dhists_1<- graph_features_to_histograms(props_1)
#'  dhists_2<- graph_features_to_histograms(props_2)
#'  # Obtain the minimum NetEMD_edges between the histograms 
#'  netemd_single_pair(dhists_1[[1]],dhists_2[[1]],method = "optimise",smoothing_window_width = 0)
#' @export
netemd_single_pair <- function(dhist1, dhist2, method = "optimise",
                               smoothing_window_width = 0) {
  # Present dhists as smoothed or unsmoothed histograms depending on the value
  # of smoothing_window_width
  # NOTE: This MUST be done prior to any variance normalisation as the
  # calculation of variance differs depending on whether or not the histograms
  # are smoothed (i.e. we need to ensure that the smoothing_window_width
  # attribute of the dhists is set to the smoothing_window_width parameter
  # provided by the caller)
  # TODO: Consider moving the smoothing of histograms outside to the user's
  # calling code. It feels a bit untidy in here.
  if (smoothing_window_width == 0) {
    dhist1 <- as_unsmoothed_dhist(dhist1)
    dhist2 <- as_unsmoothed_dhist(dhist2)
  } else {
    dhist1 <- as_smoothed_dhist(dhist1, smoothing_window_width)
    dhist2 <- as_smoothed_dhist(dhist2, smoothing_window_width)
  }
  
  # Store means and variances to calculate offset later
  mean1 <- dhist_mean_location(dhist1)
  mean2 <- dhist_mean_location(dhist2)
  
  var1 <- dhist_variance(dhist1)
  var2 <- dhist_variance(dhist2)
  
  # Mean centre histograms. This helps with numerical stability as, after
  # variance normalisation, the differences between locations are often small.
  # We want to avoid calculating small differences between large numbers as
  # floating point precision issues can result in accumulating inaccuracies.
  # Mean-centering histograms results in variance normalised locations being
  # clustered around zero, rather than some potentially large mean location.
  dhist1 <- mean_centre_dhist(dhist1)
  dhist2 <- mean_centre_dhist(dhist2)
  
  # Normalise histogram to unit mass and unit variance
  dhist1_norm <- normalise_dhist_variance(normalise_dhist_mass(dhist1))
  dhist2_norm <- normalise_dhist_variance(normalise_dhist_mass(dhist2))
  
  # Calculate minimal EMD
  result <- min_emd(dhist1_norm, dhist2_norm, method = method)
  # As we mean-centred the histograms prior to passing to min_emd(), the offset
  # returned is not the "true" offset for the supplied histograms. We report
  # this as the "standardised" offset.
  result$min_offset_std <- result$min_offset
  # We report the "true" offset as the offset with no mean-centring, so need to
  # adjust to reverse the earlier mean-centring
  result$min_offset <- result$min_offset + mean2 - mean1
  return(result)
}
