# HISTOGRAM FUNCTIONS
#' Discrete histogram constructor
#' 
#' Creates a discrete histogram object of class \code{dhist}, with bin 
#' \code{locations} and \code{masses} set to the 1D numeric vectors provided.
#' @param locations A 1D numeric vector specifying the discrete locations
#' of the histogram bins
#' @param masses A 1D numeric vector specifying the mass present at each 
#' location
#' @param smoothing_window_width If greater than 0, the discrete histogram will
#' be treated as having the mass at each location "smoothed" uniformly across
#' a bin centred on the location and having width = \code{smoothing_window_width}
#' (default = \code{0} - no smoothing)
#' @param sorted Whether or not to return a discrete histogram with locations 
#' and masses sorted by ascending mass (default = \code{TRUE})
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
dhist <- function(locations, masses, smoothing_window_width = 0, sorted = TRUE) {
  if(!is_numeric_vector_1d(locations)) {
    stop("Bin locations must be provided as a 1D numeric vector")
  }
  if(!is_numeric_vector_1d(masses)) {
    stop("Bin masses must be provided as a 1D numeric vector")
  }
  if(length(locations) != length(masses)) {
    stop("The number of bin locations and masses provided must be equal")
  }
  dhist <- list(locations = locations, masses = masses, 
                smoothing_window_width = smoothing_window_width)
  class(dhist) <- "dhist"
  if(sorted == TRUE) {
    dhist <- sort_dhist(dhist)
  }
  return(dhist)
}

#' Compare dhists
#' 
#' Compares all fields of the dhist and only returns treu if they are all the
#' same in both dhists
#' @param dhist1 A discrete histogram as a \code{dhist} object
#' @param dhist2 A discrete histogram as a \code{dhist} object
`==.dhist` <- function(dhist1, dhist2) {
  class(dhist1) == class(dhist2) &&
    all(mapply(`==`, dhist1$locations, dhist2$locations)) &&
    all(mapply(`==`, dhist1$masses, dhist2$masses)) && 
    dhist1$smoothing_window_width == dhist2$smoothing_window_width
}

update_dhist <- 
  function(dhist, locations = dhist$locations, masses = dhist$masses,
           smoothing_window_width = dhist$smoothing_window_width) {
    dhist$locations <- locations
    dhist$masses <- masses
    dhist$smoothing_window_width <- smoothing_window_width
    return(dhist)
    }

#' Set dhist smoothing
#' 
#' Returns a "smoothed" copy of a \code{dhist} object with its 
#' \code{smoothing_window_width} attribute set to the value provided 
#' \code{smoothing_window_width} parameter.
#' @param dhist A discrete histogram as a \code{dhist} object
#' @param smoothing_window_width If greater than 0, the discrete histogram will
#' be treated as having the mass at each location "smoothed" uniformly across
#' a bin centred on the location and having width = \code{smoothing_window_width}
#' @return A copy of a \code{dhist} object with its \code{smoothing_window_width}
#' attribute set  to the value provided \code{smoothing_window_width} parameter.
#' @export
as_smoothed_dhist <- function(dhist, smoothing_window_width) {
  dhist <- update_dhist(dhist, smoothing_window_width = smoothing_window_width)
  return(dhist)
}

#' Remove dhist smoothing
#' 
#' Returns an "unsmoothed" copy of a \code{dhist} object with its 
#' \code{smoothing_window_width} attribute set to 0.
#' @param dhist A discrete histogram as a \code{dhist} object
#' @return A copy of a \code{dhist} object with its \code{smoothing_window_width}
#' attribute set to 0.
#' @export
as_unsmoothed_dhist <- function(dhist) {
  dhist <- update_dhist(dhist, smoothing_window_width = 0)
  return(dhist)
}

#' Check if an object is a \code{dhist} discrete histogram
#' 
#' Checks if the input object is of class \code{dhist}. If \code{fast_check} is 
#' \code{TRUE} then the only check is whether the object has a class attribute of 
#' \code{dhist}. If \code{fast_check} is \code{FALSE} (default), then checks
#' are also made to ensure that the object has the structure required of a
#' \code{dhist} object. 
#' @param x An arbitrary object
#' @param fast_check Boolean flag indicating whether to perform only a 
#' superficial fast check limited to checking the object's class attribute 
#' is set to \code{dhist} (default = \code{TRUE})
#' @export
is_dhist <- function(x, fast_check = TRUE) {
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
dhist_from_obsSLOW <- function(observations) {
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
    if (any(is.na(observations))) {
        stop("NA observed in features")
    }
  results <- counts_from_observations(matrix(observations))
  # Construct histogram object
  hist <- dhist(locations = results[,1], masses = results[,2])
  return(hist)
}

#' Generate interpolating empirical cumulative mass function (ECMF) for 
#' a discrete histogram
#' 
#' @param dhist A discrete histogram as a \code{dhist} object
#' @return An interpolating ECMF as an \code{approxfun} object. This function
#' will return the interpolated cumulative mass for a vector of arbitrary 
#' locations. If \code{dhist$smoothing_window_width} is zero, the ECMF will be
#' piecewise constant. If \code{dhist$smoothing_window_width} is one, the ECMF
#' will be piece-wise linear. If \code{dhist$smoothing_window_width} is any
#' other value, the ECMF will not correctly represent the cumulative mass.
#' @export
dhist_ecmf <- function(dhist) {
  # Ensure histogram is sorted in order of increasing location
  dhist <- sort_dhist(dhist, decreasing = FALSE)
  # Determine cumulative mass at each location
  cum_mass <- cumsum(dhist$masses)
  # Generate ECMF
  if(dhist$smoothing_window_width == 0) {
    # Avoid any issues with floating point equality comparison completely when
    # no smoothing is occurring
    x_knots <- dhist$locations
    interpolation_method <- "constant"
  } else {
    hw <- dhist$smoothing_window_width / 2
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
    cum_mass_lower <- c(0, utils::head(cum_mass_upper, num_locs -1))
    # 3. Identify upper limits within machine precision of the lower limit of 
    # the next location
    diff <- abs(utils::head(upper_limits, num_locs -1) - 
                  utils::tail(lower_limits, num_locs -1))
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
    # where gradient of ECMF may change (i.e. the "knots" of the ECMF)
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
  dhist_ecmf <- stats::approxfun(x = x_knots, y = cum_mass, 
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
#' @return x_knots A list of "knots" for the ECMF, containing all x-values at 
#' which the y-value changes gradient (i.e. the x-values between which the ECMF
#' does its constant or linear interpolation)
#' @export
ecmf_knots <- function(dhist_ecmf) {
  # dhist_ecmf is a stats::approxfun object and is either a piecewise constant
  # or piece-wise linear function, with the x argument of the underlying 
  # approxfun set to the inflexion points (or knots) of the pricewise function
  # We simply recover the value of the x argument by evaluating "x" in the 
  # environment of the dhist_ecmf approxfun
  eval(expression(x), envir=environment(dhist_ecmf))
}

#' Calculate area between two discrete histogram empirical cumulative 
#' mass functions (ECMFs)
#' 
#' @param dhist_ecmf1 An object of class \code{dhist_ecmf}, returned from a call 
#' to the \code{dhist_ecmf} function
#' @param dhist_ecmf2 An object of class \code{dhist_ecmf}, returned from a call 
#' to the \code{dhist_ecmf} function
#' @return area The area between the two discrete histogram ECMFs, calculated as
#' the integral of the absolute difference between the two ECMFs
#' @export
area_between_dhist_ecmfs <- function(dhist_ecmf1, dhist_ecmf2) {
  # Ensure ECMFs have compatible types
  ecmf_type1 <- attr(dhist_ecmf1, "type")
  ecmf_type2 <- attr(dhist_ecmf2, "type")
  if(ecmf_type1 != ecmf_type2) {
    stop("ECMFs must have the same type")
  }
  ecmf_type <- ecmf_type1
  # Determine all possible locations where either ECMF changes gradient ("knots")
  x1 <- ecmf_knots(dhist_ecmf1)
  x2 <- ecmf_knots(dhist_ecmf2)
  x <- sort(union(x1, x2))
  # Calculate the cumulative density at each of these locations for both ECMFs
  ecm1 <- dhist_ecmf1(x)
  ecm2 <- dhist_ecmf2(x)
  # Set some other values used in either case
  num_segs <- length(x) - 1
  x_lower <- utils::head(x, num_segs)
  x_upper <- utils::tail(x, num_segs)
  # Depending on the ECDF type, we calculate the area between ECMFs differently
  if(ecmf_type == "constant") {
    # Area of each rectangular segment between ECMFs is the absolute difference
    # between the ECMFs at the lower limit of the segment * the width of the
    # segement
    ecm_diff <- abs(ecm2 - ecm1)
    ecm_diff_lower <- utils::head(ecm_diff, num_segs)
    segment_width <- abs(x_upper - x_lower)
    segment_areas <- ecm_diff_lower * segment_width
  } else if(ecmf_type == "linear") {
    # --------------------------------------------------------------
    # Determine area between pairs of linear segments from each ECMF
    # --------------------------------------------------------------
    y1_lower <- utils::head(ecm1, num_segs)
    y1_upper <- utils::tail(ecm1, num_segs)
    y2_lower <- utils::head(ecm2, num_segs)
    y2_upper <- utils::tail(ecm2, num_segs)
    # Determine if ECMFs intersect within each segment. The linear segments from
    # each ECMF will only intersect if the ordering of the y-components of their
    # start and end endpoints are different (i.e. the ECMF with the y-component
    # at the start of the segment has the higher y-component at the end of the 
    # segment). An equivalent expression of this condition is that the signs 
    # of the differences between the y-components of the two linear ECMF
    # segments will differ at the start (lower x-bound) and end (upper x-bound)
    # of a segment
    y_diff_lower <- y2_lower - y1_lower
    y_diff_upper <- y2_upper - y1_upper
    # To check for opposing signs, we check if the product of the y-component
    # differences is less than 0. This formulation classifies cases where either
    # one (triangle) or both (coincident) differences are 0 as trapeziums for
    # the purposes of area calculation. The fast bowtie area calculation we use
    # will try and divide by zero if the ECMF segments are coincident, while the
    # trapezium area formula will give the correct area both for triangles and
    # coincident lines.
    bowtie <- (y_diff_lower * y_diff_upper) < 0
    trapezium <- !bowtie
    # Set segment areas
    x_diff <- x_upper - x_lower
    segment_areas <- rep(NaN, num_segs)
    # Use bowtie area for bowties
    segment_areas[bowtie] <- 
      segment_area_bowtie(x_diff = x_diff[bowtie], 
                          y_diff_lower = y_diff_lower[bowtie],
                          y_diff_upper = y_diff_upper[bowtie])
    # Use trapezium area for other shapes (trapeziums, triangles and zero-area
    # co-linear)
    segment_areas[trapezium] <- 
      segment_area_trapezium(x_diff = x_diff[trapezium],
                             y_diff_lower = y_diff_lower[trapezium],
                             y_diff_upper = y_diff_upper[trapezium])
  } else {
    stop("ECMF type not recognised")
  }
  area <- sum(segment_areas)
  return(area)
}

segment_area_trapezium <- function(x_diff, y_diff_lower, y_diff_upper) {
  height_trapezium <- abs(x_diff)
  base_trapezium <- abs(y_diff_lower)
  top_trapezium <- abs(y_diff_upper)
  segment_area <- 0.5 * height_trapezium * (base_trapezium + top_trapezium)
}

segment_area_bowtie <- function(x_diff, y_diff_lower, y_diff_upper) {
  # Note that this formula only holds when y_diff_lower and y_diff_upper have
  # opposite signs and are not both zero.
  # See issue #21 for verification that this approach is equivalent to the
  # previous approach when the above conditions hold.
  segment_area <- 0.5 * x_diff * (y_diff_lower^2 + y_diff_upper^2) / 
    (abs(y_diff_lower) + abs(y_diff_upper))
}

#' Area between two offset Empirical Cumulative Mass Functions (ECMFs)
#' 
#' @param ecmf1 An Empirical Cululative Mass Function (ECMF) object of class 
#' \code{dhist_ecmf}
#' @param ecmf2 An Empirical Cululative Mass Function (ECMF) object of class 
#' \code{dhist_ecmf}
#' @param offset An offset to add to all locations of the first ECMF. Postive
#' offsets will shift the ECMF to the right and negative ones to the left.
#' @return area The area between the two ECMFs, calculated as the integral of 
#' the absolute difference between the two ECMFs
area_between_offset_ecmfs <- function(ecmf1, ecmf2, offset) {
  # Construct ECMFs for each normalised histogram
  ecmf1 <- dhist_ecmf(shift_dhist(dhist1_norm, offset))
  ecmf2 <- dhist_ecmf(dhist2_norm)
  area_between_dhist_ecmfs(ecmf1, ecmf2)
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
  mean_centred_locations <- dhist$locations - dhist_mean_location(dhist)
  # Variance is E[X^2] - E[X]. However, for mean-centred data, E[X] is zero, 
  # so variance is simply E[X^2]. Centring prior to squaring also helps avoid
  # any potential integer overfloww issues (R uses a signed 32-bit integer 
  # representation, so cannot handle integers over ~2.1 billion)
  if(dhist$smoothing_window_width == 0) {
    # For unsmoothed discrete histograms, the mass associated with each location
    # is located precisely at the lcoation. Therefore cariance (i.e. E[X^2])
    # is the mass-weighted sum of the mean-centred locations
    variance <- sum(dhist$masses * (mean_centred_locations)^2) / sum(dhist$masses)
  } else {
    # For smoothed histograms, the mass associated with each location is "smoothed"
    # uniformly across a bin centred on the location with width = smoothing_window_width
    # Variance (i.e. E[X^2]) is therefore the mass-weighted sum of the integrals
    # of x^2 over the mean-centred bins at each location.
    hw = dhist$smoothing_window_width / 2
    bin_lowers <- mean_centred_locations - hw
    bin_uppers <- mean_centred_locations + hw
    # See comment in issue #21 on Github repository for verification that E[X^2]
    # is calculated as below for a uniform bin
    bin_x2_integrals <- (bin_lowers^2 + bin_uppers^2 + bin_lowers*bin_uppers) / 3
    variance <- sum(dhist$masses * bin_x2_integrals) / sum(dhist$masses)
  }
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
  dhist <- update_dhist(dhist,locations = centred_locations)
  return(dhist)
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
  dhist <- update_dhist(dhist, masses = normalised_masses)
  return(dhist)
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
  # Special case for histograms with only one location and no smoothing. 
  # Variance is zero / undefined so normalisation fails. Just return bin centres
  # unchanged
  if(length(dhist$locations) == 1 && dhist$smoothing_window_width == 0) {
    dhist <- dhist
  } else {
    # Centre locations on mean, divide centred locations by standard deviation
    # then uncentre them
    std_dev <- dhist_std(dhist)
    centred_locations <- (dhist$locations - dhist_mean_location(dhist))
    normalised_centred_locations <- centred_locations / std_dev
    normalised_locations <- normalised_centred_locations + dhist_mean_location(dhist)
    dhist <- update_dhist(dhist, locations = normalised_locations)
    # If smoothing_window_width not zero, then update it to reflect the variance
    # normalisation
    if(dhist$smoothing_window_width != 0) {
      normalised_smoothing_window_width <- dhist$smoothing_window_width / std_dev
      dhist <- update_dhist(dhist, smoothing_window_width = normalised_smoothing_window_width)
    }
  }
  return(dhist)
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
  dhist1 <- update_dhist(dhist1, locations = locations1, masses = masses1)
  dhist2 <- update_dhist(dhist2, locations = locations2, masses = masses2)
  return(list(dhist1 = dhist1, dhist2 = dhist2))
}

#' Check if 1D numeric vector
#' 
#' Check if a variable is a 1D numeric vector by checking that:
#' \itemize{
#'   \item \code{is_numeric(input)}: Input is vector, matrix, array or list of numbers
#'   \item \code{is_null(dim(input))}: Input is not a matrix or array
#' }
#' @param input Arbitrary object
#' @return TRUE if input is a 1D numeric vector. FALSE otherwise.
is_numeric_vector_1d <- function(input) {
  return(is.numeric(input) && purrr::is_null(dim(input)))
}
