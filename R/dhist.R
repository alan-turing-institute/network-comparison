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
  dhist <- sort_dhist(dhist)
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
    # Early return if fast check requested
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

#' Generate interpolating empirical cumulative mass function (ECMF) for 
#' a discrete histogram
#' 
#' @param dhist A discrete histogram as a \code{dhist} object
#' @param smoothing_window_width Width of "top-hat" smoothing window to apply to
#' "smear" point masses across a finite width in the real domain prior to 
#' calculating the ECMF. This will result in a piecewise linear ECMF rather than
#' the stepwise constant ECMF obtained with no smoothing. Default is 0, which 
#' results in no smoothing and a stepwise constant ECMF. Care should be taken 
#' to select a \code{smoothing_window_width} that is appropriate for the 
#' discrete domain (e.g.for the integer domain a width of 1 is the natural
#' choice)
#' @return An interpolating ECMF as an \code{approxfun} object. This function
#' will return the interpolated cumulative mass for a vector of arbitrary locations.
#' @export
dhist_ecmf <- function(dhist, smoothing_window_width = 0, normalise_mass = FALSE, normalise_variance = FALSE) {
  # Normalise histogram to unit mass if requested
  if(normalise_mass) {
    dhist <- normalise_dhist_mass(dhist)
  }
  # Normalise histogram to unit variance if requested
  if(normalise_variance) {
    dhist <- normalise_dhist_variance(dhist)
  }
  # Ensure histogram is sorted in order of increasing location
  dhist <- sort_dhist(dhist, decreasing = FALSE)
  # Determine cumulative mass at each location
  cum_mass <- cumsum(dhist$masses)
  # Generate ECDM
  if(smoothing_window_width == 0) {
    # Avoid any issues with floating point equality comparison completely when
    # no smoothing is occurring
    x_knots <- dhist$locations
    interpolation_method <- "constant"
  } else {
    hw <- smoothing_window_width / 2
    # Determine set of "knots" that define the ECMF and the value of the ECMF
    # at each knot
    # 1. Initial knot candidates are at +/- half the smoothing window width
    # from the discrete locations
    num_locs <- length(dhist$locations)
    lower_limits <- dhist$locations - hw
    upper_limits <- dhist$locations + hw
    cum_mass_lower <- cum_mass
    cum_mass_upper <- cum_mass
    # 2. Set lower limit cumulative masses to have the same value as the 
    # upper limit of the previous location. This ensures constant interpolation
    # between the upper limit of one location and the lower limit of the next
    cum_mass_lower <- c(0, head(cum_mass_upper, num_locs -1))
    # 3. Identify upper limits within machine precision of the lower limit of 
    # the next location
    diff <- abs(head(upper_limits, num_locs -1) - tail(lower_limits, num_locs -1))
    tolerance <- .Machine$double.eps
    drop_indexes <- which(diff <= tolerance)
    # 4. Drop upper limits and associated cumulative masses where a lower 
    # limit exists at the same location (to within machine precision).
    # NOTE: We need to skip this step entirely if there are no upper limits to
    # drop as vector[-0] returns an empty vector rather than all entries in the
    # vector.
    if(length(drop_indexes) > 0) {
      upper_limits <- upper_limits[-drop_indexes]
      cum_mass_upper <- cum_mass_upper[-drop_indexes]
    }
    # 5. Combine location limits and cumulative masses to define all points
    # where gradient of ECDM may change (i.e. the "knots" of the ECMF)
    x_knots <- c(lower_limits, upper_limits)
    cum_mass <- c(cum_mass_lower, cum_mass_upper)
    # 6. Sort knots and cumulative mass by increasing order of location
    sorted_indexes <- sort(x_knots, decreasing = FALSE, index.return = TRUE)$ix
    x_knots <- x_knots[sorted_indexes]
    cum_mass <- cum_mass[sorted_indexes]
    # 7. Set interpolation method
    interpolation_method <- "linear"
  }
  # Construct ECMF
  max_mass <- cum_mass[length(cum_mass)]
  dhist_ecmf <- approxfun(x = x_knots, y = cum_mass, 
                          method = interpolation_method, yleft = 0, 
                          yright = max_mass, f = 0, ties = min)
  class(dhist_ecmf) <- c("dhist_ecmf", class(dhist_ecmf))
  attr(dhist_ecmf, "type") <- interpolation_method
  return(dhist_ecmf)
}

#' Get "knots" for discrete histogram empirical cumulative mass function
#' (ECMF). The "knots" are the x-values at which the y-value of the ECDM changes
#' gradient (i.e. the x-values between which the ECMF does its constant or 
#' linear interpolates)
#' 
#' @param dhist_ecmf An object of class \code{dhist_ecmf}, returned from a call 
#' to the \code{dhist_ecmf} function
#' @param smoothing_window_width Width of "top-hat" smoothing window to apply to
#' @return x_knots A list of "knots" for the ECMF, containing all x-values at 
#' which the y-value changes gradient (i.e. the x-values between which the ECMF
#' does its constant or linear interpolation)
#' @export
knots.dhist_ecmf <- function(dhist_ecmf, ...) {
  eval(expression(x), envir=environment(dhist_ecmf))
}
#' Sort discrete histogram
#' 
#' Sort a discrete histogram so that locations are in increasing (default) or 
#' decreasing order
#' @param dhist A discrete histogram as a \code{dhist} object
#' @param decreasing Logical indicating whether histograms should be sorted in 
#' increasing (default) or decreasing order of location
#' @export
sort_dhist <- function(dhist, decreasing = FALSE) {
  sorted_indexes <- sort(dhist$locations, decreasing = decreasing, index.return = TRUE)$ix
  dhist$masses <- dhist$masses[sorted_indexes]
  dhist$locations <- dhist$locations[sorted_indexes]
  return(dhist)
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
#' NOTE: Does not apply bias correction (i.e. N-1 denominator) as masses 
#' may not represent counts so N is not necessarily known
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
#' NOTE: Does not apply bias correction (i.e. N-1 denominator) as masses 
#' may not represent counts so N is not necessarily known
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
#' @export
normalise_dhist_mass <- function(dhist) {
  total_mass <- sum(dhist$masses)
  normalised_masses <- dhist$masses / total_mass
  return(dhist(masses = normalised_masses, locations = dhist$locations))
}

#' Normalise a discrete histogram to unit variance
#' 
#' Normalises a discrete histogram to unit variance by dividing each centred
#' location by the standard deviation of the discrete histogram before 
#' decentering
#' @param dhist A discrete histogram as a \code{dhist} object
#' @return A discrete histogram normalised to have variance 1
#' @export
normalise_dhist_variance <- function(dhist) {
  # Special case for histograms with only one location. Variance is zero / undefined
  # so normalisation fails. Just return bin centres unchanged
  if(length(dhist$locations) == 1) {
    return(dhist)
  }
  centred_locations <- (dhist$locations - dhist_mean_location(dhist))
  normalised_centred_locations <- centred_locations / dhist_std(dhist)
  normalised_locations <- normalised_centred_locations + dhist_mean_location(dhist)
  return(dhist(masses = dhist$masses, locations = normalised_locations))
}

#' Harmonise a pair of discrete histograms to share a common set of locations
#' 
#' Where a location only exists in one histogram, add this location to the other
#' histogram with zero mass. This ensures that all location exist in both 
#' histograms.
#' @param dhist1 A discrete histogram as a \code{dhist} object
#' @param dhist2 A discrete histogram as a \code{dhist} object
#' @return Harmonised histograms
#' @export
harmonise_dhist_locations <- function(dhist1, dhist2) {
  # Identify missing locations in each histogram
  missing_locations1 <- setdiff(dhist2$locations, dhist1$locations)
  missing_locations2 <- setdiff(dhist1$locations, dhist2$locations)
  # Add missing locations to end of each histogram
  locations1 <- c(dhist1$locations, missing_locations1)
  locations2 <- c(dhist2$locations, missing_locations2)
  # Assign these extra locations zero mass
  masses1 <- c(dhist1$masses, rep(0, length(missing_locations1)))
  masses2 <- c(dhist2$masses, rep(0, length(missing_locations2)))
  # Construct a new histogram using the dhist constructor to ensure that the
  # harmonised histograms have the same properties as if they had been 
  # constructed with the additional bins in the first place
  # (e.g. sorted by location)
  dhist1 <- dhist(locations = locations1, masses = masses1)
  dhist2 <- dhist(locations = locations2, masses = masses2)
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
