#' NetEMDs between all graph pairs using provided Graphlet-based Degree 
#' Distributions 
#' @param gdds List containing sets of Graphlet-based Degree Distributions for 
#' all graphs being compared
#' @param method The method to use to find the minimum EMD across all potential 
#' offsets for each pair of histograms. Default is "optimise" to use
#' R's built-in \code{optimise} method to efficiently find the offset with the 
#' minimal EMD. However, this is not guaranteed to find the global minimum if 
#' multiple local minima EMDs exist. You can alternatively specify the 
#' "fixed_step" method, which will exhaustively evaluate the EMD between the 
#' histograms at overlapping offsets separated by a fixed step. The size of the 
#' fixed step is 1/2 the the minimum spacing between locations in either
#' histogram after normalising to unit variance
#' @param step_size Additional optional argument for "fixed_step" method allowing
#' user to specify their own minumum step size. Note that this step size is 
#' applied to the histograms after they have been normalised to unit variance
#' @param return_details Logical indicating whether to return the individual
#' minimal EMDs and associated offsets for all pairs of histograms
#' @param smoothing_window_width Width of "top-hat" smoothing window to apply to
#' "smear" point masses across a finite width in the real domain. Default is 0, 
#' which  results in no smoothing. Care should be taken to select a 
#' \code{smoothing_window_width} that is appropriate for the discrete domain 
#' (e.g.for the integer domain a width of 1 is the natural choice)
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
net_emds_for_all_graphs <- function(
  gdds, method = "optimise", step_size = NULL, smoothing_window_width = 0, 
  return_details = FALSE, mc.cores = getOption("mc.cores", 2L)) {
  comp_spec <- graph_cross_comparison_spec(gdds)
  out <- purrr::simplify(parallel::mcmapply(function(index_a, index_b) {net_emd(
    gdds[[index_a]], gdds[[index_b]], method = method, return_details = return_details,
    smoothing_window_width = smoothing_window_width)
    }, comp_spec$index_a, comp_spec$index_b, SIMPLIFY = FALSE))
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

graph_cross_comparison_spec <- function(gdds) {
  indexes <- as.data.frame(t(utils::combn(1:length(gdds),2)))
  names <- as.data.frame(cbind(names(gdds)[indexes[,1]], names(gdds)[indexes[,2]]))
  spec <- cbind(names, indexes)
  colnames(spec) <- c("name_a", "name_b", "index_a", "index_b")
  return(spec)
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
#' R's built-in \code{optimise} method to efficiently find the offset with the 
#' minimal EMD. However, this is not guaranteed to find the global minimum if 
#' multiple local minima EMDs exist. You can alternatively specify the 
#' "fixed_step" method, which will exhaustively evaluate the EMD between the 
#' histograms at overlapping offsets separated by a fixed step. The size of the 
#' fixed step is 1/2 the the minimum spacing between locations in either
#' histogram after normalising to unit variance
#' @param step_size Additional optional argument for "fixed_step" method allowing
#' user to specify their own minumum step size. Note that this step size is 
#' applied to the histograms after they have been normalised to unit variance
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
net_emd <- function(dhists1, dhists2, method = "optimise", step_size = NULL, 
                    return_details = FALSE, smoothing_window_width = 0) {
  # Require either a pair of "dhist" discrete histograms or two lists of "dhist"
  # discrete histograms
  pair_of_dhist_lists <- all(purrr::map_lgl(dhists1, is_dhist)) && all(purrr::map_lgl(dhists2, is_dhist))
  
  # If input is two lists of "dhist" discrete histograms, determine the minimum
  # EMD and associated offset for pairs of histograms taken from the same 
  # position in each list
  if(pair_of_dhist_lists) {
    details <- purrr::map2(dhists1, dhists2, function(dhist1, dhist2) {
      net_emd_single_pair(dhist1, dhist2, method = method, step_size = step_size,
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
    return(net_emd(list(dhists1), list(dhists2), method = method, 
                   step_size = step_size, return_details = return_details,
                   smoothing_window_width = smoothing_window_width))
  }
}

net_emd_single_pair <- function(dhist1, dhist2, method = "optimise", 
                                step_size = NULL, smoothing_window_width = 0) {
  # Require input to be a pair of "dhist" discrete histograms 
  if(!(is_dhist(dhist1) && is_dhist(dhist2))) {
    stop("All inputs must be 'dhist' discrete histogram objects")
  }
  # Set default method if not supplied
  if(missing(method)) {
    method = "optimise"
  }
  
  # Present dhists as smoothed or unsmoothed histograms depending on the value 
  # of smoothing_window_width
  # NOTE: This MUST be done prior to any variance normalisation as the 
  # calculation of variance differes depending on whether or not the histograms 
  # are smoothed (i.e. we need to ensure that the smoothing_window_width 
  # attribute of the dhists is set to the smoothing_window_width parameter
  # provided by the caller)
  # TODO: Consider moving the smoothing of histograms outside to the user's 
  # calling code. It feels a bit untidy in here.
  if(smoothing_window_width == 0) {
    dhist1 <- as_unsmoothed_dhist(dhist1)
    dhist2 <- as_unsmoothed_dhist(dhist2)
  } else {
    dhist1 <- as_smoothed_dhist(dhist1, smoothing_window_width)
    dhist2 <- as_smoothed_dhist(dhist2, smoothing_window_width)
  }
  
  # Normalise histogram to unit mass and unit variance
  dhist1_norm <- normalise_dhist_variance(normalise_dhist_mass(dhist1))
  dhist2_norm <- normalise_dhist_variance(normalise_dhist_mass(dhist2))
  
  # Determine minimum and maximum offset of range in which histograms overlap
  # (based on sliding histogram 1)
  min_offset <- min(dhist2_norm$locations) - max(dhist1_norm$locations)
  max_offset <- max(dhist2_norm$locations) - min(dhist1_norm$locations)

  if(method == "fixed_step" && is.null(step_size)) {
    # Set default step_size for "fixed_step" method if step_size not provided
    location_spacing <- function(l) {
      l <- sort(l)
      tail(l, length(l)-1) - head(l, length(l)-1)
    }
    min_location_sep1 <- min(location_spacing(dhist1_norm$locations))
    min_location_sep2 <- min(location_spacing(dhist2_norm$locations))
    step_size <- min(min_location_sep1, min_location_sep2)/2
  }
  
  emd_offset <- function(offset) {
    # Construct ECMFs for each normalised histogram
    ecmf1 <- dhist_ecmf(shift_dhist(dhist1_norm, offset))
    ecmf2 <- dhist_ecmf(dhist2_norm)
    area_between_dhist_ecmfs(ecmf1, ecmf2)
  }
  
  # Define optimise method for picking minimal EMD offset
  # 1. "optimise" method
  min_emd_opt <- function() {
    # Set lower and upper range for optimise algorithm to be somewhat wider than
    # range defined by the minimum and maximum offset. This guards against a
    # couple of issues that arise if the optimise range is exactly min_offset 
    # to max_offset
    # 1) If lower and upper are equal, the optimise method will throw and error
    # 2) It seems that optimise is not guaranteed to explore its lower and upper
    #    bounds, even in the case where one of them is the offset with minimum
    #    EMD
    buffer <- 0.1
    soln <- optimise(emd_offset, lower = (min_offset -buffer), upper = (max_offset + buffer))
    min_emd <- soln$objective
    min_offset <- soln$minimum
    return(list(min_emd = min_emd, min_offset = min_offset))
  }
  # 2. "fixed_step" method
  min_emd_step <- function(step_size) {
    offsets <- seq(min_offset, max_offset, by = step_size)
    emds <- purrr::map_dbl(offsets, emd_offset)
    min_idx <- which.min(emds)
    min_emd <- emds[min_idx]
    min_offset <- offsets[min_idx]
    return(list(min_emd = min_emd, min_offset = min_offset))
  }
   
  # Determine minimum EMD across all offsets
  min_emd_details <- switch(EXPR = method, 
         optimise = min_emd_opt(),
         fixed_step = min_emd_step(step_size),
         stop("Supplied 'method' not recognised")
  )
  return(min_emd_details)
}

#' Earth Mover's Distance (EMD) 
#' 
#' Calculates the Earth Mover's Distance (EMD) between two discrete histograms
#' @param dhist1 A \code{dhist} discrete histogram object
#' @param dhist2 A \code{dhist} discrete histogram object
#' @return Earth Mover's Distance between the two discrete histograms
#' @export
emd <- function(dhist1, dhist2) {
  # Require the inputs to be "dhist" objects
  if(!(is_dhist(dhist1) && is_dhist(dhist2))) {
    stop("All inputs must be 'dhist' discrete histogram objects")
  }
  # Use efficient difference of cumulative histogram method that can also 
  # handle non-integer bin masses and location differences
  emd_cs(dhist1, dhist2)
}

#' Earth Mover's Distance (EMD) using linear programming (LP)
#' 
#' Takes two sets of histogram bin masses and bin centres and calculates the 
#' Earth Mover's Distance between the two histograms by solving the Transport
#' Problem using linear programming.
#' 
#' WARNING: Linear Programming solution will only give a correct answer if all
#' masses and distances between bin centres are integers.
#' @param bin_masses1 Bin masses for histogram 1
#' @param bin_masses2 Bin masses for histogram 2
#' @param bin_centres1 Bin centres for histogram 1
#' @param bin_centres2 Bin centres for histogram 2
#' @return Earth Mover's Distance between the two input histograms
#' @export
emd_lp <- function(bin_masses1, bin_masses2, bin_centres1, bin_centres2) {
  num_bins1 <- length(bin_masses1)
  num_bins2 <- length(bin_masses2)
  
  # Check inputs: All bins in each histogram must have a mass and centre, so
  # the bin_mass and bin_centre vectors for each histogram must have the same
  # length.
  if(length(bin_centres1) != num_bins1) {
    stop("Number of bin masses and bin centres provided for histogram 1 must be equal")
  }
  if(length(bin_centres2) != num_bins2) {
    stop("Number of bin masses and bin centres provided for histogram 2 must be equal")
  }
  
  # Generate cost matrix
  cost_mat <- cost_matrix(bin_centres1, bin_centres2)
 
  # Linear Programming solver requires all bin masses and transportation costs 
  # to be integers to generate correct answer
  if(!isTRUE(all.equal(bin_masses1, floor(bin_masses1)))) {
    stop("All bin masses for histogram 1 must be integers for accurate Linear Programming solution")
  }
  if(!isTRUE(all.equal(bin_masses2, floor(bin_masses2)))) {
    stop("All bin masses for histogram 2 must be integers for accurate Linear Programming solution")
  }
  if(!isTRUE(all.equal(cost_mat, floor(cost_mat)))) {
    stop("All costs must be integers for accurate Linear Programming solution")
  } 
  row_signs <- rep("==", num_bins1)
  col_signs <- rep("<=", num_bins2)
  s <- lpSolve::lp.transport(cost.mat = cost_mat, row.signs = row_signs, 
                    col.signs = col_signs, row.rhs = bin_masses1, 
                    col.rhs = bin_masses2)
  return(s$objval)
}

#' Earth Mover's Distance (EMD) using the difference of cumulative sums method
#' 
#' Takes two discrete histograms and calculates the Wasserstein / Earth Mover's
#' Distance between the two histograms by summing the absolute difference 
#' between the two cumulative histograms.
#' @references 
#' Calculation of the Wasserstein Distance Between Probability Distributions on the Line
#' S. S. Vallender, Theory of Probability & Its Applications 1974 18:4, 784-786
#' \url{http://dx.doi.org/10.1137/1118101}
#' @param dhist1 A discrete histogram as a \code{dhist} object
#' @param dhist2 A discrete histogram as a \code{dhist} object
#' @return Earth Mover's Distance between the two input histograms
#' @export
emd_cs <- function(dhist1, dhist2) {
  # Generate Empirical Cumulative Mass Functions (ECMFs) for each discrete histogram
  ecmf1 <- dhist_ecmf(dhist1)
  ecmf2 <- dhist_ecmf(dhist2)
  # Calculate the area between the two ECMFs
  area <- area_between_dhist_ecmfs(ecmf1, ecmf2)
  return(area)
}

#' Inter-bin cost matrix from bin centres
#' 
#' Generates a matrix for the cost of moving a unit of mass between each bin in 
#' histogram 1 and each bin in histogram 2.
#' @param bin_centres1 Bin centres for histogram 1
#' @param bin_centres2 Bin centres for histogram 2
#' @return Cost matrix
cost_matrix <- function(bin_centres1, bin_centres2) {
  # Calculate distances between all bins in network 1 and all bins in network 2
  num_bins1 <- length(bin_centres1) 
  num_bins2 <- length(bin_centres2)
  loc_mat1 <- matrix(bin_centres1, nrow = num_bins1, ncol = num_bins2, byrow = FALSE)
  loc_mat2 <- matrix(bin_centres2, nrow = num_bins1, ncol = num_bins2, byrow = TRUE)
  cost_mat <- abs(loc_mat1 - loc_mat2)
  return(cost_mat)
}
