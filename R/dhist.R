

# HISTOGRAM FUNCTIONS
#' Discrete histogram consructor
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
    stop("The number of bin locations and masses provided must be equals")
  }
  dhist <- list(locations = locations, masses = masses, class = "hist")
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

#' Check if 1D numeric vector
#' 
#' Check if a variable is a 1D numeric vector by checking that:
#' \itemize{
#'   \item \code{is_numeric(input)}: Input is vector, matrix, array or list of numbers
#'   \item \code{is_null(dim(input))}: Input is not a matrix or array
#' }
#' @return TRUE if input is a 1D numeric vector. FALSE otherwise.
#' @export
is_numeric_vector_1d <- function(input) {
  return(is_numeric(input) && is_null(dim(input)))
}
