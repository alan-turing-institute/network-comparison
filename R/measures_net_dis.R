#' @export
count_graphlets_ego_scaled <- function(graph, max_graphlet_size, 
                                            neighbourhood_size) {
  # Calculate ego-network graphlet counts
  ego_graphlet_counts <- 
    count_graphlets_ego(graph, max_graphlet_size = max_graphlet_size, 
                        neighbourhood_size = neighbourhood_size)
  # Calculate total number of k-tuples in ego-network for each graphlet type
  ego_networks <- igraph::make_ego_graph(graph, order = neighbourhood_size)
  ego_node_counts <- purrr::map_int(ego_networks, function(g) {length(igraph::V(g))})
  graphlet_key <- graphlet_key(max_graphlet_size)
  ego_ktuples <- simplify2array(purrr::map(graphlet_key$node_count,
     function(k) { choose(ego_node_counts, k) }))
  # Scale ego-network graphlet counts by dividing by total number of k-tuples in
  # ego-network (where k is graphlet size)
  ego_graphlet_counts_scaled <- ego_graphlet_counts / ego_ktuples
  return(ego_graphlet_counts_scaled)
}

#' @export
interval_indexes <- function(x, breaks) {
  # Get indexes for the intervals each value of x falls into. Setting 
  # all.inside = TRUE ensures that the minimum and maximum values of x will be 
  # assigned to the intervals they bound.
  indexes <- findInterval(x, breaks, all.inside = TRUE)
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
  # There is one less interval than there are breaks
  num_intervals <- length(breaks) - 1
  # Get indexes for the intervals each value of x falls into.
  x_interval_indexes <- interval_indexes(x, breaks)
  # Find the lowest interval with fewer than the minimum required count.
  # Not all intervals are guaranteed to have members in x. If they don't, they 
  # won't appear in x_interval_indexes. We therefore append the full list of 
  # indexes prior to counting and subtract 1 from all counts afterwards to get 
  # an accurate count that includes indexes with no members with zero counts
  all_interval_indexes <- 1:num_intervals
  interval_index_counts <- plyr::count(c(x_interval_indexes, all_interval_indexes))
  interval_index_counts$freq <- interval_index_counts$freq - 1
  
  # Find the first interval with fewer members than the minimum specified count
  merge_position <- Position(function(i) i < min_count, interval_index_counts$freq)
  # Not all intervals are guaranteed to have members, so convert the index 
  # provided by Position into an index into the full interval list and then add
  merge_interval_index <- interval_index_counts$x[merge_position]
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