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
dhist_ecmf <- function(dhist, smoothing_window_width = 0) {
  # Ensure histogram is sorted in order of increasing location
  dhist <- sort_dhist(dhist, decreasing = FALSE)
  # Determine cumulative mass at each location
  cum_mass <- cumsum(dhist$masses)
  # Generate ECMF
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
  x1 <- knots(dhist_ecmf1)
  x2 <- knots(dhist_ecmf2)
  x <- sort(union(x1, x2))
  # Calculate the cumulative density at each of these locations for both ECMFs
  ecm1 <- dhist_ecmf1(x)
  ecm2 <- dhist_ecmf2(x)
  # Calculate the spacing between x-values and the difference between ECMFs at
  # each x-value (used in all cases)
  ecm_diff <- abs(ecm2 - ecm1)
  num_segs <- length(ecm_diff) - 1
  x_lower <- head(x, num_segs)
  x_upper <- tail(x, num_segs)
  segment_width <- abs(x_upper - x_lower)
  # Depending on the ECDF type, we calculate the area between ECMFs differently
  if(ecmf_type == "constant") {
    # Area of each rectangular segment between ECMFs is the absolute difference
    # between the ECMFs at the lower limit of the segment * the width of the
    # segement
    ecm_diff_lower <- head(ecm_diff, num_segs)
    segment_areas <- ecm_diff_lower * segment_width
  } else if(ecmf_type == "linear") {
    # --------------------------------------------------------------
    # Determine area between pairs of linear segments from each ECMF
    # --------------------------------------------------------------
    x_lowers <- head(x, length(x) - 1)
    x_uppers <- tail(x, length(x) - 1)
    segment_areas <- purrr::map2_dbl(x_lowers, x_uppers, function(x_l, x_u) {
      segment_area_piecewise_linear(x_l = x_l, x_u = x_u, f1 = dhist_ecmf1, f2 = dhist_ecmf2)
    })
  } else {
    stop("ECMF type not recognised")
  }
  area <- sum(segment_areas)
  return(area)
}

segment_area_piecewise_linear <- function(f1, f2, x_l, x_u) {
  y1_l <- f1(x_l)
  y1_u <- f1(x_u)
  line1 <- line(point(x_l, y1_l), point(x_u, y1_u))
  y2_l <- f2(x_l)
  y2_u <- f2(x_u)
  line2 <- line(point(x_l, y2_l), point(x_u, y2_u))
  intersection <- segment_intersection(line1, line2)
  if(intersection$type == "intersecting") {
    # Segments form two triangles in a "bow-tie" (or potentially a single
    # triangle if segments are just touching, but this can be covered with 
    # the same area formula)
    height_lower_triangle <- intersection$point$x - x_l
    height_upper_triangle <- x_u - intersection$point$x
    base_lower_triangle <- abs(y2_l - y1_l)
    base_upper_triangle <- abs(y2_u - y1_u)
    area_lower_triangle <- 0.5 * base_lower_triangle * height_lower_triangle
    area_upper_triangle <- 0.5 * base_upper_triangle * height_upper_triangle
    segment_area <- area_lower_triangle + area_upper_triangle
  } else if(intersection$type == "parallel" 
            || intersection$type == "non-intersecting") {
    # Segments 
    top_trapezium <- abs(y2_l - y1_l)
    base_trapezium <- abs(y2_u - y1_u)
    height_trapezium <- abs(x_u - x_l)
    area_trapezium <- 0.5 * (top_trapezium + base_trapezium) * height_trapezium
    segment_area <- area_trapezium
  } else {
    stop("Invalid intersection type")
  }
  return(segment_area)
}

point <- function(x, y) {
  list(x = x, y = y)
}

line <- function(point1, point2) {
  list(point1 = point1, point2 = point2)
}

segment_intersection <- function(line1, line2) {
  # Implementation of line segment intersection algorithm adaptedfrom 
  # "Computational Geometry in C", J. O'Rourke, 1994, pp 220-226
  # http://crtl-i.com/PDF/comp_c.pdf
  
  # Translate our line segment endpoints into the representation used by O'Rourke
  a <- line1$point1; a0 <- a$x; a1 <- a$y
  b <- line1$point2; b0 <- b$x; b1 <- b$y
  c <- line2$point1; c0 <- c$x; c1 <- c$y
  d <- line2$point2; d0 <- d$x; d1 <- d$y
  
  # Generate s, t and D
  # Generate s*D and t*D, so we can do our within segment check prior to division
  # NOTE: There is a missing '-' sign on p 221 that should multiply the definition
  # of 't' in equation 7.2 by -1 (as is seen in the sample code). This has been 
  # added here
  sD <- (a0 * (d1 - c1) + c0 * (a1 - d1) + d0 * (c1 - a1))
  tD <- -(a0 * (c1 - b1) + b0 * (a1 - c1) + c0 * (b1 - a1))
  D <- a0 * (d1 - c1) + b0 * (c1 - d1) + d0 * (b1 - a1) + c0 * (a1 - b1)
  s <- sD / D
  t <- tD / D
  
  if(D == 0) {
    # Lines are parallel
    type <- "parallel"
    point <- point(NaN, NaN)
  }  else {
    # We make comparisons between sD, tD and D, rather than between s, t and 1
    # as in O'Rourke. This is to avoid dividing by potentially small numbers 
    # prior to these comparisons. However, this means we need to consider the
    # cases when sD, tD are both positive and negative.
    sD_positive_between_0_and_D <- ((0 <= sD) && (sD <= D))
    sD_negative_between_0_and_D <- ((0 >= sD) && (sD >= D))
    tD_positive_between_0_and_D <- ((0 <= tD) && (tD <= D))
    tD_negative_between_0_and_D <- ((0 >= tD) && (tD >= D))
    if((sD_positive_between_0_and_D || sD_negative_between_0_and_D ) &&
       (tD_positive_between_0_and_D || tD_negative_between_0_and_D)) {
      # The values of sD and tD will lie within the segment, so the lines
      # will intersect within the segment
      type <- "intersecting"
    } else {
      type <- "non-intersecting"
    }
    x <- a0 + s * (b0 - a0)
    y <- a1 + s * (b1 - a1)
    point <- point(x, y)
  }
  return(list(type = type, point = point))
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
