
#' Generate graph cross-comparison specification
#' 
#' @param graph_edges List of graphs in ORCA-compatible graph edge format
#' @return M
#' @export
graph_cross_comparison_spec <- function(graph_edges) {
  indexes <- as.data.frame(t(utils::combn(1:length(graph_edges),2)))
  names <- as.data.frame(cbind(attr(graph_edges, "name")[indexes[,1]], attr(graph_edges, "name")[indexes[,2]]))
  spec <- cbind(names, indexes)
  colnames(spec) <- c("name_a", "name_b", "index_a", "index_b")
  return(spec)
}

# WIP function to cross-comexpand.grid(i = 1:3, j = 4:6) pare a set of edge lists
# cross_compare_graphs <- function(graph_edges, feature_spec = "orbits1") {
#   spec <- graph_cross_comparison_spec(graph_edges)
# }
