# GRAPH FUNCTIONS
library(igraph)

#' Integer index edge list from igraph
#' 
#' Takes a igraph graph object and generates an edgelist where each edge is
#' represented by the integer IDs of its vertices. Note that, where a graph
#' has isolated vertices, the IDs for these vertices will not be present
#' in the edge list. Where a graph has no isolated vertices, the edge list will 
#' include all vertex IDs from 1 to numVertices.
#' @param graph igraph graph object
#' @return edges A 2 x numEdges edgelist with integer vertex indices
#' @export
integer_index_edge_list <- function(graph) {
  # Use igraph method to get edge list with edges specified using vertex ID 
  # (indexes) rather than names
  edges <- igraph::get.edgelist(graph, names = FALSE)
  # Convert edge list from numeric to integer
  edges <- structure(vapply(edges, as.integer, integer(1)), dim = dim(edges))
  return(edges)
}

#' Make arbitrary igraph ORCA compatible
#' 
#' Takes a igraph graph object and does the following to make the graph 
#' compatible with the ORCA fast orbit calculation package:
#'   1. Makes the graph undirected
#'   2. Removes loops (where both endpoints of an edge are the same vertex)
#'   3. Removes multiple edges (i.e. ensuring only one edge exists for each 
#'      pair of endpoints)
#'   4. Removes isolated vertices (i.e. vertices with no edges after the 
#'      previous alterations)
#' @param graph Original igraph graph object
#' @return graph Processed igraph graph object
#' @export
orca_graph <- function(graph) {
  # Ensure graph is undirected
  graph <- igraph::as.undirected(graph)
  # Remove loops (where both endpoints of an edge are the same vertex) and 
  # multiple edges (where two edges have the same endpoints [in the same order
  # for directed graphs])
  graph <- igraph::simplify(graph, remove.loops = TRUE, remove.multiple = TRUE)
  # Remove vertices that have no edges to ensure no missing IDs when generating
  # an ID-based edge list from the graph
  isolated_vertex_indices <- (igraph::degree(graph) == 0)
  graph <- igraph::delete.vertices(graph, isolated_vertex_indices)
  return(graph)
}

#' Load graph from file and process to be ORCA compatible
#' 
#' Loads a graph from file into an igraph graph object and does the following 
#' to make the graph compatible with the ORCA fast orbit calculation package:
#'   1. Makes the graph undirected
#'   2. Removes loops (where both endpoints of an edge are the same vertex)
#'   3. Removes multiple edges (i.e. ensuring only one edge exists for each 
#'      pair of endpoints)
#'   4. Removes isolated vertices (i.e. vertices with no edges after the 
#'      previous alterations)
#' @param file Path to graph file
#' @param format Format of graph file
#' @return graph ORCA compatible igraph graph object
#' @export
read_orca_graph <- function(file, format = "ncol") {
  # Read graph from file
  graph <- igraph::read.graph(file = file, format = format)
  # Make graph compatibe with ORCA fast orbit calculation package
  graph <- orca_graph(graph)
  return(graph)
}
