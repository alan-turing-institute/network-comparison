library("lpSolve")
library("purrr")

#' NetEMD Network Earth Mover's Distance
#' 
#' Calculates the minimum Earth Mover's Distance (EMD) between two histograms
#' after normalising to unit mass and variance.
#' This is calculated as follows:
#'   1. Normalise each histogram to have unit mass and unit variance
#'   2. "Slide" histogram 1 over histogram 2 by varying the bin offset between 
#'      the two histograms, calculating the EMD at each offset.
#'   3. The NetEMD is the minimum EMD observed across all offsets
#' @param bin_masses1 Bin masses for histogram 1
#' @param bin_masses2 Bin masses for histogram 2
#' @param bin_centres1 Bin centres for histogram 1
#' @param bin_centres2 Bin centres for histogram 2
#' @return NetEMD
#' @export
net_emd <- function(bin_masses1, bin_masses2, bin_centres1, bin_centres2, step = 1) {
  # Normalise histograms to unit variance
  bin_centres1 <- normalise_histogram_variance(bin_masses1, bin_centres1)
  bin_centres2 <- normalise_histogram_variance(bin_masses2, bin_centres2)
  # Normalise histograms to unit mass
  bin_masses1 <- normalise_histogram_mass(bin_masses1)
  bin_masses2 <- normalise_histogram_mass(bin_masses2)
  
  # Determine minimum and maximum offset of range in which histograms overlap
  # if sliding histogram 1
  min_offset <- min(bin_centres2) - max(bin_centres1)
  max_offset <- max(bin_centres2) - min(bin_centres1)

  # Offset histogram 1 by all possible overlapping offsets and calculate EMD at
  # each offset, keeping track of minimum EMD across all offsets. This is the 
  # NetEMD
  offsets <- seq(min_offset, max_offset, by = step)

  emds <- purrr:::map_dbl(offsets, function(offset) {emd_cs(bin_masses1, bin_masses2, bin_centres1 + offset, bin_centres2)})
  net_emd <- min(emds)
  return(net_emd)
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
 
  # Linear Programming solv%er requires all bin masses and transportation costs 
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
#' Takes two sets of histogram bin masses and bin centres and calculates the 
#' Earth Mover's Distance between the two histograms by summing the absolute 
#' difference between the two cumulative histograms.
#' @references 
#' Calculation of the Wasserstein Distance Between Probability Distributions on the Line
#' S. S. Vallender, Theory of Probability & Its Applications 1974 18:4, 784-786
#' \url{http://dx.doi.org/10.1137/1118101}
#' @param bin_masses1 Bin masses for histogram 1
#' @param bin_masses2 Bin masses for histogram 2
#' @param bin_centres1 Bin centres for histogram 1
#' @param bin_centres2 Bin centres for histogram 2
#' @return Earth Mover's Distance between the two input histograms
#' @export
emd_cs <- function(bin_masses1, bin_masses2, bin_centres1, bin_centres2) {
  
  # Check inputs: All bins in each histogram must have a mass and centre, so
  # the bin_mass and bin_centre vectors for each histogram must have the same
  # length.
  num_bins1 <- length(bin_masses1)
  num_bins2 <- length(bin_masses2)
  if(length(bin_centres1) != num_bins1) {
    stop("Number of bin masses and bin centres provided for histogram 1 must be equal")
  }
  if(length(bin_centres2) != num_bins2) {
    stop("Number of bin masses and bin centres provided for histogram 2 must be equal")
  }
  
  # Ensure both histograms have entries for all bins that appear in either histogram
  ah <- augment_histograms(bin_masses1, bin_masses2, bin_centres1, bin_centres2)
  bin_masses1 <- ah$bin_masses1
  bin_masses2 <- ah$bin_masses2
  bin_centres1 <- ah$bin_centres1
  bin_centres2 <- ah$bin_centres2
  
  # Sort bin centres and masses in bin centre order
  sorted_indexes1 <- sort(bin_centres1, decreasing = FALSE, index.return = TRUE)$ix
  sorted_indexes2 <- sort(bin_centres2, decreasing = FALSE, index.return = TRUE)$ix
  sorted_masses1 <- bin_masses1[sorted_indexes1]
  sorted_masses2 <- bin_masses2[sorted_indexes2]
  sorted_centres1 <- bin_centres1[sorted_indexes1]
  sorted_centres2 <- bin_centres2[sorted_indexes2]
  
  # Generate cumulative histogram masses
  cum_mass1 <- cumsum(sorted_masses1)
  cum_mass2 <- cumsum(sorted_masses2)
  
  # Determine spacing between bin centres
  vector_spacing <- function(v) {
    tail(v, length(v)-1) - head(v, length(v)-1)
  }
  sorted_centre_spacing <- vector_spacing(sorted_centres1)
  
  # Discrete integration of absolute difference between cumulative histograms
  # - Get difference between cumulative histograms at each bin with data
  cum_mass_diff_at_bins <- apply(rbind(cum_mass1, cum_mass2), 2, function(x) abs(x[1] - x[2]))
  # - Multiply the difference at each bin by the distance to the next bin to get
  # - the area under the difference curve between consecutive bins (bins are 
  # - not guaranteed to be equally sized or equally spaced and bins with no 
  # - mass in either histogram may not be present at all)
  cum_mass_diff_area_between_bins <- 
    head(cum_mass_diff_at_bins, length(cum_mass_diff_at_bins) - 1) * 
    sorted_centre_spacing
  
  return(sum(cum_mass_diff_area_between_bins))
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

#' Augment a pair of histograms to share a common set of bins
#' 
#' Where a bin centre only exists in one histogram, add a zero mass bin with the
#' same centre to the other histogram to ensure that all bins exist in both 
#' histograms. This makes some hsitogram operations easier.
#' @param bin_masses1 Bin masses for histogram 1
#' @param bin_masses2 Bin masses for histogram 2
#' @param bin_centres1 Bin centres for histogram 1
#' @param bin_centres2 Bin centres for histogram 2
#' @return Augmented histograms
augment_histograms <- function(bin_masses1, bin_masses2, bin_centres1, bin_centres2) {
  # Identify missing bin centres in each histogram
  missing_centres1 <- setdiff(bin_centres2, bin_centres1)
  missing_centres2 <- setdiff(bin_centres1, bin_centres2)
  # Add missing centres to end of each histogram
  bin_centres1 <- c(bin_centres1, missing_centres1)
  bin_centres2 <- c(bin_centres2, missing_centres2)
  # Assign these extra bins zero mass
  bin_masses1 <- c(bin_masses1, rep(0, length(missing_centres1)))
  bin_masses2 <- c(bin_masses2, rep(0, length(missing_centres2)))
  
  list(bin_masses1 = bin_masses1, bin_masses2 = bin_masses2, 
       bin_centres1 = bin_centres1, bin_centres2 = bin_centres2)
}
#' Normalise histogram to unit mass
#' 
#' Normalises histogram to unit mass by dividing each bin mass by the total of 
#' the non-normalised bin masses
#' @param bin_masses Bin masses for histogram
#' @return Bin masses normalised to sum to 1
normalise_histogram_mass <- function(bin_masses) {
  total_mass <- sum(bin_masses)
  normalised_masses <- bin_masses / total_mass
  return(normalised_masses)
}

#' Normalise histogram to unit variance
#' 
#' Normalises histogram to unit variance by dividing each bin centre by the 
#' standard deviation of the hsitogram
#' @param bin_masses Bin masses for histogram
#' @param bin_centres Bin centres for histogram
#' @return Bin centres normalised to give a histogram of variance 1
normalise_histogram_variance <- function(bin_masses, bin_centres) {
  mean_centre <- sum(bin_masses * bin_centres) / sum(bin_masses)
  centred_centres <- (bin_centres - mean_centre)
  normalised_centred_centres <- centred_centres / histogram_std(bin_masses, bin_centres)
  normalised_bin_centres <- normalised_centred_centres + mean_centre
  return(normalised_bin_centres)
}

#' Calculate histogram variance
#' 
#' Calculates variance directly from the histogram by using bin centres weighted
#' by bin masses. Where a histogram represents a quantised summary of an 
#' underlying sample, it is more accurate to calculate variance directly from 
#' the sample data if available.
#' NOTE: Does not apply bias correction (i.e. N-1 denominator) as bin_masses 
#' may not represent bin counts so N is not necessarily known
#' @param bin_masses Bin masses for histogram
#' @param bin_centres Bin centres for histogram
#' @return Variance of histogram
histogram_variance <- function(bin_masses, bin_centres) {
  mean_centre <- sum(bin_masses * bin_centres) / sum(bin_masses)
  variance <- sum(bin_masses * (bin_centres - mean_centre)^2) / sum(bin_masses)
  return(variance)
}

#' Calculate histogram standard deviation
#' 
#' Calculates standard deviation directly from the histogram by using bin 
#' centres weighted by bin masses. Where a histogram represents a quantised 
#' summary of an underlying sample, it is more accurate to calculate standard
#' deviation directly from the sample data if available. 
#' NOTE: Does not apply bias correction (i.e. N-1 denominator) as bin_masses 
#' may not represent bin counts so N is not necessarily known
#' @param bin_masses Bin masses for histogram
#' @param bin_centres Bin centres for histogram
#' @return Standard deviation of histogram
histogram_std <- function(bin_masses, bin_centres) {
  return(sqrt(histogram_variance(bin_masses, bin_centres)))
}
