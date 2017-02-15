library("purrr")

# HISTOGRAM FUNCTIONS
#' Discrete histogram constructor
#' 
#' Creates a discrete histogram object of class \code{dhist}, with bin 
#' \code{locations} and \code{masses} set to the 1D numeric vectors provided.
#' @param \code{locations} A 1D numeric vector specifying the discrete locations
#' of the histogram bins
#' @param \code{masses} A 1D numeric vector specifying the mass present at each 
#' location
#' @return A sparse discrete histogram. Format is a \code{dhist} object, which
#' is a list of class \code{dhist} with the following named elements:
#' \itemize{
#'   \item \code{locations}: A 1D numeric vector of discrete locations
#'   \item \code{masses}: A 1D numeric vector of the mass present at each location
#' }
#' Note that locations where no mass is present are not included in the returned
#' \code{dhist} object. Mass in these discrete histograms is treated as being
#' present precisely at the specified location. Discrete histograms should not be used
#' for data where observations have been grouped into bins representing ranges 
#' of observation values.
#' @export
dhist <- function(locations, masses) {
  if(!is_numeric_vector_1d(locations)) {
    stop("Bin locations must be provided as a 1D numeric vector")
  }
  if(!is_numeric_vector_1d(masses)) {
    stop("Bin masses must be provided as a 1D numeric vector")
  }
  if(length(locations) != length(masses)) {
    stop("The number of bin locations and masses provided must be equal")
  }
  dhist <- list(locations = locations, masses = masses)
  class(dhist) <- "dhist"
  return(dhist)
}

#' Check if an object is a \code{dhist} discrete histogram
#' 
#' Checks if the input object is of class \code{dhist}. If \code{fast_check} is 
#' \code{TRUE} then the only check is whether the object has a class attribute of 
#' \code{dhist}. If \code{fast_check} is \code{FALSE} (default), then checks
#' are also made to ensure that the object has the structure required of a
#' \code{dhist} object. 
#' @param \code{x} An arbitrary object
#' @param \code{fast_check} Boolean flag inficating whether to perform only a 
#' superficial fast check limited to checking the object's class attribute 
#' is set to \code{dhist} (default = \code{FALSE})
#' @export
is_dhist <- function(x, fast_check = FALSE) {
  # Quick check that relies on user not to construct variables with "dhist" class
  # that do not have the required elements
  has_class_attr <-(class(x) == "dhist")
  if(fast_check) {
    # Early return is fast check requested
    return(has_class_attr)
  }
  # Otherwise check structure
  has_locations <- purrr::contains(attr(x, "name"), "locations")
  has_masses <- purrr::contains(attr(x, "name"), "masses")
  # Require list with correct class and presence of 1D numeric vector named 
  # elements "locations" and "masses"
  return(has_class_attr
         && purrr::is_list(x)
         && has_locations
         && has_masses
         && is_numeric_vector_1d(x$locations)
         && is_numeric_vector_1d(x$masses))
}

#' Discrete histogram from observations
#' 
#' Generate a sparse discrete histogram from a set of discrete numeric observations
#' @param observations A vector of discrete numeric observations
#' @return A sparse discrete histogram. Format is a \code{dhist} object, which
#' is a list of class \code{dhist} with the following named elements:
#' \itemize{
#'   \item \code{locations}: A 1D numeric vector of discrete locations
#'   \item \code{masses}: A 1D numeric vector of the mass present at each location
#' }
#' @export
dhist_from_obs <- function(observations) {
  # Require 1D numeric vector
  if(!is_numeric_vector_1d(observations)) {
    stop("Observations must be provided as a 1D numeric vector")
  }

  # Identify unique observations
  locations <- sort(unique(observations))
  # Count occurences of each unique obervation
  counts <- sapply(locations, function(location) {sum(observations == location)})
  # Construct histogram object
  hist <- dhist(locations = locations, masses = counts)
  return(hist)
}

#' Shift discrete histogram
#' 
#' Shift the locations of a discrete histogram rightwards on the x-axis by the 
#' specified amount
#' @param dhist A discrete histogram as a \code{dhist} object
#' @param shift The distance to add to all locations
#' @return A shifted discrete histogram as a \code{dhist} object
#' @export
shift_dhist <- function(dhist, shift) {
  dhist$locations <- (dhist$locations + shift)
  return(dhist)
}

#' Calculate mean location for a discrete histogram
#' 
#' Calculates mean location for a discrete histogram by taking a weighted sum
#' of each location weighted by the fraction of the total histogram mass at that
#' location.
#' @param dhist A discrete histogram as a \code{dhist} object
#' @return The mass-weighted mean location
#' @export
dhist_mean_location <- function(dhist) {
  sum((dhist$masses/ sum(dhist$masses)) * dhist$locations) 
}

#' Calculate variance of a discrete histogram
#' 
#' Calculates variance directly from the discrete histogram by using locations
#' weighted by  masses. 
#' NOTE: Does not apply bias correction (i.e. N-1 denominator) as bin_masses 
#' may not represent bin counts so N is not necessarily known
#' @param dhist A discrete histogram as a \code{dhist} object
#' @return Variance of histogram
#' @export
dhist_variance <- function(dhist) {
  variance <- sum(dhist$masses * (dhist$locations - dhist_mean_location(dhist))^2) / sum(dhist$masses)
  return(variance)
}

#' Calculate standard deviation of a discrete histogram
#' 
#' Calculates standard deviation directly from the discrete histogram by using 
#' locations weighted by masses.
#' NOTE: Does not apply bias correction (i.e. N-1 denominator) as bin_masses 
#' may not represent bin counts so N is not necessarily known
#' @param dhist A discrete histogram as a \code{dhist} object
#' @return Standard deviation of histogram
#' @export
dhist_std <- function(dhist) {
  return(sqrt(dhist_variance(dhist)))
}

#' Centre a discrete histogram around its mean location
#' 
#' Centres a discrete histogram around its mass-weighted mean location by 
#' subtracting the mass-weighted mean from each location.
#' @param dhist A discrete histogram as a \code{dhist} object
#' @return The mass-weighted mean location
#' @export
mean_centre_dhist <- function(dhist) {
  centred_locations <- dhist$locations - dhist_mean_location(dhist)
  return(dhist(masses = dhist$masses, locations = centred_locations))
}

#' Normalise a discrete histogram to unit mass
#' 
#' Normalises a discrete histogram to unit mass by dividing each mass by the 
#' total of the non-normalised masses
#' @param dhist A discrete histogram as a \code{dhist} object
#' @return A discrete histogram normalised to have mass 1
normalise_dhist_mass <- function(dhist) {
  total_mass <- sum(dhist$masses)
  normalised_masses <- dhist$masses / total_mass
  return(dhist(masses = normalised_masses, locations = dhist$locations))
}

#' Normalise a discrete histogram to unit variance
#' 
#' Normalises a discrete histogram to unit variance by dividing each location by
#' the standard deviation of the discrete histogram
#' @param dhist A discrete histogram as a \code{dhist} object
#' @return A discrete histogram normalised to have variance 1
normalise_dhist_variance <- function(dhist) {
  # Special case for histograms with only one bin. Variance is zero / undefined
  # so normalisation fails. Just return bin centres unchanged
  if(length(dhist$locations) == 1) {
    return(dhist)
  }
  centred_locations <- (dhist$locations - dhist_mean_location(dhist))
  normalised_centred_locations <- centred_locations / dhist_std(dhist)
  normalised_locations <- normalised_centred_locations + dhist_mean_location(dhist)
  return(dhist(masses = dhist$masses, locations = normalised_locations))
}

#' Homogonise a pair of discrete histograms to share a common set of locations
#' 
#' Where a location only exists in one histogram, add this location to the other
#' histogram with zero mass. This ensures that all location exist in both 
#' histograms.
#' @param dhist1 A discrete histogram as a \code{dhist} object
#' @param dhist2 A discrete histogram as a \code{dhist} object
#' @return Augmented histograms
harmonise_dhist_locations <- function(dhist1, dhist2) {
  # Identify missing locations in each histogram
  missing_locations1 <- setdiff(dhist2$locations, dhist1$locations)
  missing_locations2 <- setdiff(dhist1$locations, dhist2$locations)
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
