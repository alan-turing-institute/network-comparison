library("purrr")

# VECTOR FUNCTIONS
rotl_vec <- function(vec, lshift) {
  num_els <- length(vec)
  select_mask <- ((1:num_els + lshift) %% num_els) 
  select_mask[select_mask==0] <- num_els
  return(vec[select_mask])
}

rotr_vec <- function(vec, rshift) {
  return(rotl_vec(vec, -rshift))
}

# HISTOGRAM FUNCTIONS
#' Discrete histogram
#' 
#' Generate a sparse discrete histogram from a set of integer observations
#' @param observations A vector of integer-valued observations
#' @return A sparse discrete histogram. Format is a list with the following 
#' named elements:
#' \itemize{
#'   \item \code{$locations}: Integer-valued locations where one or more observations 
#'   exist
#'   \item \code{$masses}: Integer counts of observations at each location
#' }
#' @export
discrete_hist <- function(observations) {
  # Require all observations to be integers
  if(!isTRUE(all.equal(observations, floor(observations)))) {
    stop("All observations must be integers to generate a discrete histogram")
  }
  
  locations <- sort(unique(observations))
  counts <- sapply(locations, function(location) {sum(observations == location)})
  
  hist <- list(locations = locations, masses = counts)
  return(hist)
}

#' Shift histogram
#' 
#' Shift the bin locations of a histogram rightwards on the x-axis by the 
#' specificed shift
#' @param histogram A histogram as a list with named members `masses` and `locations`
#' @param shift The distance to add to all bin locations
#' @return A  shifted histogram as a list with named members `masses` and `locations`
#' @export
shift_histogram <- function(histogram, shift) {
  histogram$locations <- histogram$locations + shift
  return(histogram)
}

