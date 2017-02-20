library("purrr")

# HISTOGRAM FUNCTIONS
#' Continuous histogram constructor
#' 
#' Creates a sparse continuous histogram object of class \code{chist}, with bin 
#' \code{lower_edges}, \code{upper_edges} and \code{densities} set to the 1D 
#' numeric vectors provided.
#' @param \code{lower_edges} A 1D numeric vector specifying the locations of the
#' lower edges of the histogram bins
#' @param \code{upper_edges} A 1D numeric vector specifying the locations of the
#' upper edges of the histogram bins
#' @param \code{masses} A 1D numeric vector specifying the density present at 
#' each location
#' @return A sparse continuous histogram. Format is a \code{chist} object, which
#' is a list of class \code{chist} with the following named elements:
#' \itemize{
#'   \item \code{lower_edges}: A 1D numeric vector of bin lower edge locations
#'   \item \code{upper_edges}: A 1D numeric vector of bin upper edge locations
#'   \item \code{densities}: A 1D numeric vector of the density present at each location
#' }
#' Note that locations where no density is present are not included in the returned
#' \code{chist} object.
#' @export
chist <- function(lower_edges, upper_edges, densities) {
  if(!is_numeric_vector_1d(lower_edges)) {
    stop("Bin lower edge locations must be provided as a 1D numeric vector")
  }
  if(!is_numeric_vector_1d(lupper_edges)) {
    stop("Bin upper edge locations must be provided as a 1D numeric vector")
  }
  if(!is_numeric_vector_1d(densities)) {
    stop("Bin densities must be provided as a 1D numeric vector")
  }
  if(length(lower_edges) != length(upper_edges) != length(densities)) {
    stop("The number of bin edges and densities provided must be equal")
  }
  chist <- list(lower_edges = lower_Edges, upper_edges = upper_edges, densities = densities)
  class(chist) <- "chist"
  return(chist)
}

#' Check if an object is a \code{chist} continuous histogram
#' 
#' Checks if the input object is of class \code{chist}. If \code{fast_check} is 
#' \code{TRUE} then the only check is whether the object has a class attribute of 
#' \code{chist}. If \code{fast_check} is \code{FALSE} (default), then checks
#' are also made to ensure that the object has the structure required of a
#' \code{chist} object. 
#' @param \code{x} An arbitrary object
#' @param \code{fast_check} Boolean flag inficating whether to perform only a 
#' superficial fast check limited to checking the object's class attribute 
#' is set to \code{chist} (default = \code{FALSE})
#' @export
is_chist <- function(x, fast_check = FALSE) {
  # Quick check that relies on user not to construct variables with "chist" class
  # that do not have the required elements
  has_class_attr <-(class(x) == "chist")
  if(fast_check) {
    # Early return is fast check requested
    return(has_class_attr)
  }
  # Otherwise check structure
  has_lower_edges <- purrr::contains(attr(x, "name"), "lower_edges")
  has_upper_edges <- purrr::contains(attr(x, "name"), "upper_edges")
  has_densities <- purrr::contains(attr(x, "name"), "densities")
  # Require list with correct class and presence of 1D numeric vector named 
  # elements "lower_edges", "upper_edges" and "densities"
  return(has_class_attr
         && purrr::is_list(x)
         && has_lower_edges
         && has_upper_edges
         && has_densities
         && is_numeric_vector_1d(x$lower_edges)
         && is_numeric_vector_1d(x$upper_edges)
         && is_numeric_vector_1d(x$densities))
}

#' Continuous histogram from discrete observations
#' 
#' Generate a sparse continuous histogram from a set of discrete numeric observations
#' @param observations A vector of discrete numeric observations
#' @param bin_width Width of all histogram bins. Required as the histogram will 
#' be sparse
#' @return A sparse continuous histogram. Format is a \code{chist} object, which
#' is a list of class \code{chist} with the following named elements:
#' \itemize{
#'   \item \code{lower_edges}: A 1D numeric vector of bin lower edge locations
#'   \item \code{upper_edges}: A 1D numeric vector of bin upper edge locations
#'   \item \code{densities}: A 1D numeric vector of the density present at each location
#' }
#' @export
chist_from_obs <- function(observations, bin_width) {
  # Require 1D numeric vector
  if(!is_numeric_vector_1d(observations)) {
    stop("Observations must be provided as a 1D numeric vector")
  }

  # Identify unique observations
  locations <- sort(unique(observations))
  location_spacing <- tail(locations, length(locations)-1) - head(locations, length(locations)-1)
  # If locations are closer than requested bin_width, we cannot construct an
  # accurate histogram, as we are constructing one bin per unique observation
  # rather than binning multiple unique observations together
  if(min(location_spacing) < bin_width) {
    stop("Minimum difference between observations exceeds requested bin width")
  }
  # Count occurences of each unique obervation
  counts <- sapply(locations, function(location) {sum(observations == location)})
  # Construct histogram object
  hw = bin_with / 2
  total_mass = sum(counts * bin_width)
  hist <- chist(lower_edges = locations - hw, upper_edges = locations + hw, densities = counts / total_mass)
  return(hist)
}

#' Shift continuous histogram
#' 
#' Shift the locations of a continuous histogram rightwards on the x-axis by the 
#' specified amount
#' @param chist A continuous histogram as a \code{chist} object
#' @param shift The distance to add to all locations
#' @return A shifted continuous histogram as a \code{chist} object
#' @export
shift_chist <- function(chist, shift) {
  chist$lower_edges <- (dhist$lower_edges + shift)
  chist$upper_edges <- (dhist$upper_edges + shift)
  return(dhist)
}

#' Calculate mean location for a continuous histogram
#' 
#' Calculates mean location for a continuous histogram by taking a weighted sum
#' of each bin midpoint weighted by the fraction of the total histogram mass in
#' that bin
#' @param chist A continuous histogram as a \code{chist} object
#' @return The mass-weighted mean location
#' @export
chist_mean_location <- function(chist) {
  midpoints <- chist$lower_edges + ((chist$upper_edges - chist$lower_edges) / 2)
  widths <- (chist$upper_edges - chist$lower_edges)
  masses <- chist$densities * widths
  sum((masses / sum(masses)) * midpoints) 
}

#' Calculate variance of a continuous histogram
#' 
#' Calculates variance directly from the continuous histogram by using bin
#' midpoints weighted by bin masses. 
#' NOTE: Does not apply bias correction (i.e. N-1 denominator) as bin masses 
#' may not represent bin counts so N is not necessarily known
#' @param dhist A continuous histogram as a \code{chist} object
#' @return Variance of histogram
#' @export
chist_variance <- function(chist) {
  # Use E[(x - E[x])^2] formulation rather than E[x^2] - E[x]^2 formulation for
  # numerical stability
  mean_location <- chist_mean_location(chist)
  masses <- chist$densities * widths
  centred_lower_edges <- chist$lower_edges - mean_location
  centred_upper_edges <- chist$upper_edges - mean_location
  variance <- sum(masses * (centred_lower_edges^2 + centred_upper_edges^2 + centred_lower_edges*centred_upper_edges)/3)
  return(variance)
}

#' Calculate standard deviation of a continuous histogram
#' 
#' Calculates standard continuous directly from the continuous histogram by using 
#' bin midpoints weighted by masses.
#' NOTE: Does not apply bias correction (i.e. N-1 denominator) as bin masses 
#' may not represent bin counts so N is not necessarily known
#' @param chist A discrete histogram as a \code{chist} object
#' @return Standard deviation of histogram
#' @export
chist_std <- function(chist) {
  return(sqrt(chist_variance(chist)))
}

#' Centre a continuous histogram around its mean location
#' 
#' Centres a continuous histogram around its mass-weighted mean location by 
#' subtracting the mass-weighted mean from all bin edges
#' @param dhist A continuous histogram as a \code{chist} object
#' @return The mass-weighted mean location
#' @export
mean_centre_chist <- function(chist) {
  mean_location <- chist_mean_location(chist)
  centred_lower_edges <- chist$lower_edges - mean_location
  centred_upper_edges <- chist$upper_edges - mean_location
  return(chist(lower_edges = centred_lower_edges, upper_edges = centred_upper_edges, densities = chist$densities))
}

#' Normalise a continuous histogram to unit mass
#' 
#' Normalises a continuous histogram to unit mass by dividing each bin mass by 
#' the total of the non-normalised bin masses
#' @param dhist A continuous histogram as a \code{chist} object
#' @return A continuous histogram normalised to have mass 1
normalise_chist_mass <- function(chist) {
  masses <- chist$densities * widths
  normalised_masses <- masses / sum(masses)
  normalised_densities <- normalised_masses / widths
  return(chist(lower_edges = chist$lower_edges, upper_edges = chist$upper_edges, densities = normalised_densities))
}

#' Normalise a continuous histogram to unit variance
#' 
#' Normalises a continuous histogram to unit variance by dividing each centred 
#' edge location by the standard deviation of the continuous histogram before
#' decentering
#' @param dhist A continuous histogram as a \code{chist} object
#' @return A continuous histogram normalised to have variance 1
normalise_chist_variance <- function(chist) {
  mean_location <- chist_mean_location(chist)
  centred_lower_edges <- chist$lower_edges - mean_location
  centred_upper_edges <- chist$upper_edges - mean_location
  std_dev <- chist_std(chist)
  normalised_centred_lower_edges <- centred_lower_edges / std_dev
  normalised_centred_upper_edges <- centred_upper_edges / std_dev
  normalised_lower_edges <- normalised_centred_lower_edges + mean_location
  normalised_upper_edges <- normalised_centred_upper_edges + mean_location
  return(chist(densities = chist$densities, lower_edges = normalised_lower_edges, upper_edges = normalised_upper_edges))
}

#' Homogonise a pair of continuous histograms to share a common set of bins
#' 
#' Where a bin only exists in one histogram, add this bin to the other
#' histogram with zero density This ensures that all bins exist in both 
#' histograms.
#' @param chist1 A continuous histogram as a \code{chist} object
#' @param chist2 A continuous histogram as a \code{chist} object
#' @return Augmented histograms
harmonise_chist_locations <- function(chist1, chist2) {
  # Identify missing bins in each histogram
  missing_bins <- setdiff(chist2$locations, chist1$locations)
  missing_locations2 <- setdiff(chist1$locations, chist2$locations)
  # Add missing locations to end of each histogram
  dhist1$locations <- c(dhist1$locations, missing_locations1)
  dhist2$locations <- c(dhist2$locations, missing_locations2)
  # Assign these extra locations zero mass
  dhist1$masses <- c(dhist1$masses, rep(0, length(missing_locations1)))
  dhist2$masses <- c(dhist2$masses, rep(0, length(missing_locations2)))
  
  return(list(dhist1 = dhist1, dhist2 = dhist2))
}

#' Check if 1D numeric vector
#' 
#' Check if a variable is a 1D numeric vector by checking that:
#' \itemize{
#'   \item \code{is_numeric(input)}: Input is vector, matrix, array or list of numbers
#'   \item \code{is_null(dim(input))}: Input is not a matrix or array
#' }
#' @return TRUE if input is a 1D numeric vector. FALSE otherwise.
is_numeric_vector_1d <- function(input) {
  return(purrr::is_numeric(input) && purrr::is_null(dim(input)))
}
