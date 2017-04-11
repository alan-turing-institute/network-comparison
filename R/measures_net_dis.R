
expected_graphlet_counts_ego_fn <- function(reference_graph) {
  ego_networks <- 
  ego_densities <-
  ego_graphlet_counts <-  
  ego_graphlet_counts_scaled <-
  density_bins <-
  binned_
  
}

#' Generate a set of breaks that attempt to be evenly spaces while ensuring each
#' interval has the specified minimum count
#' 
#' Starts by binning the variable by the breaks provided in \code{breaks} (if
#' \code{breaks} is a vector), or generating a set of \code{breaks} at uniformly
#' spaced intervals (if \code{breaks} is a single number). It then iteratively 
#' merges intervals with counts lower than \code{min_count} by removing breaks 
#' until all remaining intervals have counts of at least \code{min_count}.
#' 
#' @param x The variable to be binned
#' @param min_count The minimum count for each bin
#' @param breaks Either a vector containing an intital set of breaks or a single
#' number indicating how many uniformly spaced intervals to use when constructing
#' the initial set of breaks. 
#'
#' @export
adaptive_breaks <- function(x, min_count, breaks) {
  if(length(breaks) == 1) {
    # Similarly to base::cut, we interpret a single number in breaks as the
    # number of intervals required and generate these evenly spaced
    min_x <- min(x)
    max_x <- max(x)
    breaks = seq(from = min_x, to = max_x, length.out = breaks + 1)
  }
  # There is one less interval than breaks
  num_intervals <- length(breaks) - 1
  # Get indexes for the intervals each value of x falls into. Setting 
  # all.inside = TRUE ensures that the minimum and maximum values of x will be 
  # assigned to the intervals they bound.
  indexes <- findInterval(x, breaks, all.inside = TRUE)
  
  # Find the lowest interval with fewer than the minimum required count
  index_counts <- plyr::count(indexes)

  merge_position <- Position(function(i) i < min_count, index_counts$freq)
  merge_interval_index <- index_counts$x[merge_position]
  if(is.na(merge_interval_index)) {
    # If all intervals have at least the minimum count, return the breaks
    return(breaks)
  } else {
    # Remove a break to merge the low count interval with one of its neighbours
    # and recursively call this function
    if (merge_interval_index == num_intervals) {
      # If low interval is last one, we can only merge with the previous interval
      # so remove lower break for low interval
      merge_break_index <- merge_interval_index
    } else {
      # In all other cases merge low interval with next inteval by removing 
      # upper breal for low interval
      merge_break_index <- merge_interval_index + 1
    }
    return(adaptive_breaks(x, min_count, breaks[-merge_break_index]))
  }
}