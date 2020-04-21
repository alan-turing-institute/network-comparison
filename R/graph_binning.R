#' binned_densities_adaptive
#'
#' Adaptive binning function guaranteeing a minimum number of entries in each
#' bin.
#' @param densities Density values to use for binning.
#' @param min_counts_per_interval Minimum count for each bin.
#' @param num_intervals Initial number of density bins to generate.
#' TODO: Remove @export prior to publishing
#' @export
binned_densities_adaptive <- function(densities,
                                      min_counts_per_interval,
                                      num_intervals) {
  breaks <- adaptive_breaks(densities,
    min_count = min_counts_per_interval,
    breaks = num_intervals
  )
  interval_indexes <- interval_index(densities,
    breaks = breaks,
    out_of_range_intervals = FALSE
  )
  list(
    densities = densities,
    interval_indexes = interval_indexes,
    breaks = breaks
  )
}

#' Bin values into intervals based on the provided breaks
#'
#' @param x The values to be binned
#' @param breaks The boundaries between bins
#' @param out_of_range_intervals If \code{TRUE}, "out of range" values lying
#' below the first break or above the last break will be assigned to additional
#' unbounded lower and upper extrema intervals. If \code{FALSE} these "out of
#' range" values will be assigned to intervals bounded by the lowest or
#' uppermost break respectively.
#' @return A vector of bin indexes, one per value provided
#' @export
interval_index <- function(x, breaks, out_of_range_intervals = FALSE) {
  # Get indexes for the intervals each value falls into. Setting
  # all.inside = TRUE ensures that the minimum and maximum values will be
  # assigned to the intervals they bound.
  findInterval(x, breaks, all.inside = !out_of_range_intervals)
}

#' Generate a set of breaks that attempt to be evenly spaced while ensuring each
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
#' number indicating how many uniformly spaced intervals to use when
#' constructing the initial set of breaks. If a single number is provided, the
#' minumum break will be the minimum value of x and the maximum break will be
#' the maximum value of x.
#'
#' @export
adaptive_breaks <- function(x, min_count, breaks) {
  if (length(breaks) == 1) {
    # Similarly to base::cut, we interpret a single number in breaks as the
    # number of intervals required and generate these evenly spaced
    min_x <- min(x)
    max_x <- max(x)
    breaks <- seq(from = min_x, to = max_x, length.out = breaks + 1)
  }
  # There is one less interval than there are breaks
  num_intervals <- length(breaks) - 1
  # Get indexes for the intervals each value of x falls into.
  x_interval_indexes <- interval_index(x, breaks)
  # Find the lowest interval with fewer than the minimum required count.
  # Not all intervals are guaranteed to have members in x. If they don't, they
  # won't appear in x_interval_indexes. We therefore append the full list of
  # indexes prior to counting and subtract 1 from all counts afterwards to get
  # an accurate count that includes indexes with no members with zero counts
  all_interval_indexes <- 1:num_intervals
  interval_index_counts <- plyr::count(
    c(x_interval_indexes, all_interval_indexes)
  )
  interval_index_counts$freq <- interval_index_counts$freq - 1

  # Find the first interval with fewer members than the minimum specified count
  merge_position <- Position(
    function(i) i < min_count,
    interval_index_counts$freq
  )
  # Not all intervals are guaranteed to have members, so convert the index
  # provided by Position into an index into the full interval list and then add
  merge_interval_index <- interval_index_counts$x[merge_position]
  if (is.na(merge_interval_index)) {
    # If all intervals have at least the minimum count, return the breaks
    return(breaks)
  } else {
    # Remove a break to merge the low count interval with one of its neighbours
    # and recursively call this function
    if (merge_interval_index == num_intervals) {
      # If low interval is last one, we can only merge with the previous
      # interval so remove lower break for low interval
      merge_break_index <- merge_interval_index
    } else {
      # In all other cases merge low interval with next inteval by removing
      # upper breal for low interval
      merge_break_index <- merge_interval_index + 1
    }
    return(adaptive_breaks(x, min_count, breaks[-merge_break_index]))
  }
}
