#' Scaled graphlet count for ego-networks
#' 
#' Scales the ego-network count for each graphlet by dividing it by the total
#' number of possible groupings of nodes in the ego-network with the same number
#' of nodes as the graphlet.
#' @param graph A connected, undirected, simple graph as an \code{igraph} object. 
#' @param max_graphlet_size Determines the maximum size of graphlets to count. 
#' Only graphlets containing up to \code{max_graphlet_size} nodes will be counted.
#' @return ORCA-format matrix containing counts of each graphlet (columns) at 
#' each vertex in the graph (rows).
#' @export
count_graphlets_ego_scaled <- function(
  graph, max_graphlet_size, neighbourhood_size, return_ego_networks = FALSE) {
  # Calculate ego-network graphlet counts, also returning the ego networks for
  # use later in function
  ego_data <- 
    count_graphlets_ego(graph, max_graphlet_size = max_graphlet_size, 
                        neighbourhood_size = neighbourhood_size, 
                        return_ego_networks = TRUE)
  ego_graphlet_counts <- ego_data$graphlet_counts
  ego_networks <- ego_data$ego_networks
  # Scale ego-network graphlet counts by dividing by total number of k-tuples in
  # ego-network (where k is graphlet size)
  ego_graphlet_tuples <- 
    count_graphlet_tuples_ego(ego_networks, max_graphlet_size = max_graphlet_size)
  ego_graphlet_counts <- scale_graphlet_count(ego_graphlet_counts, ego_graphlet_tuples)
  # Return either graphlet counts, or graphlet counts and ego_networks
  if(return_ego_networks) {
    return(list(graphlet_counts = ego_graphlet_counts, 
                ego_networks = ego_networks))
  } else {
    return(ego_graphlet_counts)
  }
}

#' @export
mean_density_binned_graphlet_counts <- function(
  graphlet_counts, density_interval_indexes) {
  # The ego network graphlet counts are an E x G matrix with rows (E) representing
  # ego networks and columns (G) representing graphlets. We want to calculate
  # the mean count for each graphlet / density bin combination, so we will
  # use tapply to average counts for each graphlet across density bins, using
  # apply to map this operation over graphlets
  mean_density_binned_graphlet_counts <- 
    apply(graphlet_counts, MARGIN = 2, function(gc) {
      tapply(gc, INDEX = density_interval_indexes, FUN = mean)})
}

#' @export
netdis_expected_graphlet_counts <- function(
  graph, max_graphlet_size, density_breaks, density_binned_reference_counts) {
  # Look up average scaled graphlet counts for graphs of similar density
  # in the reference graph
  query_density <- igraph::edge_density(graph)
  matched_density_index <- interval_index(query_density, density_breaks)
  matched_reference_counts <- density_binned_reference_counts[matched_density_index,]
  # Scale reference counts by multiplying the reference count for each graphlet
  # by the number of possible sets of k nodes in the query graph, where k is the
  # number of nodes in the graphlet
  matched_reference_counts * count_graphlet_tuples(graph, max_graphlet_size)
}

#' @export
netdis_expected_graphlet_counts_ego <- function(
  graph, max_graphlet_size, neighbourhood_size,
  density_breaks, density_binned_reference_counts) {
  ego_networks <- make_named_ego_graph(graph, neighbourhood_size)
  expected_graphlet_counts <- purrr::map(ego_networks,
             netdis_expected_graphlet_counts,
             max_graphlet_size = max_graphlet_size,
             density_breaks = density_breaks,
             density_binned_reference_counts = density_binned_reference_counts)
  names(expected_graphlet_counts) <- names(ego_networks)
  expected_graphlet_counts
}

#' @export
zeros_to_ones <- function(v) {
  zero_index <- which(v == 0)
  v[zero_index] <- 1
  v
}

#' @export
scale_graphlet_count <- function(graphlet_count, graphlet_tuples) {
  # Avoid divide by zero errors by replacing all zeros with ones in the
  # divisor
  graphlet_count / zeros_to_ones(graphlet_tuples)
}

#' @export
count_graphlet_tuples_ego <- function(ego_networks, max_graphlet_size) {
  graphlet_tuple_counts <- 
    t(simplify2array(purrr::map(ego_networks, count_graphlet_tuples, 
             max_graphlet_size = max_graphlet_size)))
  graphlet_tuple_counts
}

#' @export
count_graphlet_tuples <- function(graph, max_graphlet_size) {
  graph_node_count <- igraph::vcount(graph)
  graphlet_key <- graphlet_key(max_graphlet_size)
  graphlet_node_counts <- graphlet_key$node_count
  graphlet_tuple_counts <- choose(graph_node_count, graphlet_node_counts)
  graphlet_tuple_counts <- setNames(graphlet_tuple_counts, graphlet_key$id)
  graphlet_tuple_counts
}

#' @export
binned_densities_adaptive <- function(densities, min_counts_per_interval, num_intervals)
{
  breaks <- adaptive_breaks(densities, min_count = min_counts_per_interval,
                            breaks = num_intervals)
  interval_indexes <- interval_index(densities, breaks = breaks,
                                     out_of_range_intervals = FALSE)
  list(densities = densities, interval_indexes = interval_indexes, breaks = breaks)
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
  x_interval_indexes <- interval_index(x, breaks)
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