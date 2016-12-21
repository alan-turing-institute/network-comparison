library("lpSolve")

#' Earth Mover's Distance (EMD) using linear programming (LP)
#' Takes two sets of feature counts and locations and calculates the 
#' Earth Mover's Distance between them. EMD is also known as the Wasserstein 
#' distance or metric
#' Required parameters are the counts of each feature for both networks
#' @param bin_counts1 Feature counts for network 1
#' @param bin_counts2 Feature counts for network 2
#' Optional parameters are the locations of each feature. If omitted, it is 
#' assumed that features are separated by unit distance
#' @param bin_locations1 Feature bin locations for network 1
#' @param bin_locations2 Feature bin locations for network 2
emd_lp <- function(bin_counts1, bin_counts2, bin_locations1, bin_locations2) {
  
  num_bins1 <- length(bin_counts1)
  if(missing(bin_locations1)) {
    bin_locations1 <- 1:num_bins1
  }
  num_locations1 <- length(bin_locations1)
  num_bins2 <- length(bin_counts2)
  if(missing(bin_locations2)) {
    bin_locations2 <- 1:num_bins2
  }
  num_locations2 <- length(bin_locations2)
  
  # Check inputs. Size of bin counts and locations must match for each histogram
  # but the histograms for the two networks can have different numbers of bins
  # as only non-zero bins are required
  if(num_bins1 != num_locations1) {
    stop("bin_counts1 and bin_locations1 must be of equal length")
  }
  if(num_bins2 != num_locations2) {
    stop("bin_counts2 and bin_locations2 must be of equal length")
  }
  
  cost_mat <- cost_matrix(bin_locations1, bin_locations2)
  row_signs <- rep("==", num_bins1)
  col_signs <- rep("<=", num_bins2)
  s <- lpSolve::lp.transport(cost.mat = cost_mat, row.signs = row_signs, 
                    col.signs = col_signs, row.rhs = bin_counts1, 
                    col.rhs = bin_counts2)
  return(s$objval)
}

#' Generates a matrix for the cost of moving a unit of mass between all pairs of
#' features. The cost is the distance between the locations of each feature pair
#' @param bin_locations1 Feature bin locations for network 1
#' @param bin_locations2 Feature bin locations for network 2 
cost_matrix <- function(bin_locations1, bin_locations2) {
  # Calculate distances between all bins in network 1 and all bins in network 2
  num_bins1 <- length(bin_locations1) 
  num_bins2 <- length(bin_locations2)
  loc_mat1 <- matrix(bin_locations1, nrow = num_bins1, ncol = num_bins2, byrow = FALSE)
  loc_mat2 <- matrix(bin_locations2, nrow = num_bins1, ncol = num_bins2, byrow = TRUE)
  cost_mat <- abs(loc_mat1 - loc_mat2)
  return(cost_mat)
}
