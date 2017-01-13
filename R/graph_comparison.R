
#' Generate graph cross-comparison specification
#' 
#' @param graph_edges List of graphs in ORCA-compatible graph edge format
#' @return M
#' @export
graph_cross_comparison_spec <- function(graph_edges) {
  names <- expand.grid(attr(graph_edges, "name"), attr(graph_edges, "name"))
  indexes <- expand.grid(1:length(graph_edges), 1:length(graph_edges))
  spec <- cbind(names, indexes)
  colnames(spec) <- c("Name A", "Name B", "Index A", "Index B")
  return(spec)
}

# WIP function to cross-compare a set of edge lists
# cross_compare_graphs <- function(graph_edges, feature_spec = "orbits1") {
#   spec <- graph_cross_comparison_spec(graph_edges)
# }
