# GRAPH FUNCTIONS
library(igraph)

#' Integer index edge list from igraph
#' 
#' Takes a igraph graph object and generates an edgelist where each edge is
#' represented by the integer indexes of its vertices. Note that, where a graph
#' has isolated vertices, the indexes for these vertices will not be present
#' in the edge list. Where a graph has no isolated vertices, the edge list will 
#' include all vertex indexes from 1 to numVertices.
#' @param graph An igraph graph object
#' @return A 2 x numEdges edgelist with vertices labelled with integer indices
#' and a 1 x numVertices node list where the Nth element of the vector contains
#' the label for the vertice represented by index N in the edgelist
#' @export
graph_to_indexed_edges <- function(graph) {
  # Use igraph method to get edge list with edges specified using vertex ID 
  # (indexes) rather than names
  edges <- igraph::get.edgelist(graph, names = FALSE)
  # Convert edge list from numeric to integer
  edges <- structure(vapply(edges, as.integer, integer(1)), dim = dim(edges))
  colnames(edges) <- c("Node A index", "Node B index")
  node_vertex_names <- igraph::get.vertex.attribute(graph, name = "name")
  attr(edges, "vertex_names") <- node_vertex_names
  return(edges)
}

#' Graph from integer index edge list
#' 
#' Takes an integer indexed edgelist (where each edge is represented by the 
#' integer indexes of its vertices) and converts it to an igraph format graph.
#' If the edge list has a "vertex_names" attribute, this will be used to name
#' the vertices in the resultant graph.
#' @param indexed_edges A 2 x numEdges edgelist with vertices labelled with 
#' integer indices, with an optional "vertex_names" attribute
#' @return An igraph graph object
#' @export
indexed_edges_to_graph <- function(indexed_edges) {
  graph <- igraph::graph_from_edgelist(indexed_edges)
  graph <- igraph::set.vertex.attribute(graph, name = "name", value = attr(indexed_edges, "vertex_names"))
  return(graph)
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
#' @return An ORCA compatible igraph graph object
#' @export
graph_to_orca_graph <- function(graph) {
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
#' @return An ORCA compatible igraph graph object
#' @export
read_orca_graph <- function(file, format = "ncol") {
  # Read graph from file
  graph <- igraph::read.graph(file = file, format = format)
  # Make graph compatibe with ORCA fast orbit calculation package
  graph <- graph_to_orca_graph(graph)
  return(graph)
}

#' Load graph from file and convert to an ORCA compatible indexed edge list
#' 
#' Loads a graph from file and converts to an indexed edge list compatible with 
#' the ORCA fast orbit calculation package by:
#'   1. Making the graph undirected
#'   2. Removing loops (where both endpoints of an edge are the same vertex)
#'   3. Removing multiple edges (i.e. ensuring only one edge exists for each 
#'      pair of endpoints)
#'   4. Removing isolated vertices (i.e. vertices with no edges after the 
#'      previous alterations)
#'   5. Relabelling vertices with an integer index from 1 to numVertices
#'   6. Converting the graph to an edge list
#' @param file Path to graph file
#' @param format Format of graph file
#' @return An ORCA compatible edge list
#' @export
read_orca_edge_list <- function(file, format = "ncol") {
  graph_to_indexed_edges(read_orca_graph(file, format))
}

#' Load all graphs in a directory, converting to ORCA compatible indexed edge lists
#' 
#' Loads graphs from all files matching the given pattern in the given directory
#' and converts them to indexed edge lists compatible with the ORCA fast orbit 
#' counting package by:
#'   1. Making the graph undirected
#'   2. Removing loops (where both endpoints of an edge are the same vertex)
#'   3. Removing multiple edges (i.e. ensuring only one edge exists for each 
#'      pair of endpoints)
#'   4. Removing isolated vertices (i.e. vertices with no edges after the 
#'      previous alterations)
#'   5. Relabelling vertices with an integer index from 1 to numVertices
#'   6. Converting the graph to an edge list
#' @param source_dir Path to graph directory
#' @param format Format of graph files
#' @param pattern Filename pattern to match graph files
#' @return A named list of ORCA compatible edge lists, with names corresponding 
#' to the names of the files each graph was loaded from
#' @export
read_all_graphs_as_orca_edge_lists <- function (source_dir, format = "ncol", pattern = ".txt") {
  # Get list of all filenames in firectory that match the pattern
  file_names <- dir(source_dir, pattern = pattern)
  # Read graph data from each ".txt" file as an ORCA-compatible indexed edge list
  edges <- purrr::map(file_names, function(file_name) {
    read_orca_edge_list(file = file.path(source_dir, file_name), format = "ncol")
  })
  # Name each edge list with the name of the file it was read from
  attr(edges, "names") <- file_names
  return(edges)
}

#' ORCA vertex graphlet orbit counts to graphlet orbit histograms
#' 
#' Converts ORCA output (counts of each graphlet orbit at each graph vertex) to 
#' a set of graphlet degree histograms (a histogram of counts across all graph 
#' vertices for each graphlet orbit) 
#' @param orca_counts ORCA output: Counts of each graphlet orbit 
#' (columns) at each graph vertex (rows)
#' @return Graphlet degree histograms: List of degree histograms for each 
#' graphlet orbit
#' @export
orca_counts_to_graphlet_orbit_degree_distribution <- function(orca_counts) {
  apply(orca_counts, 2, dhist_from_obs)
}

#' Graphlet-based Degree Distributions (GDDs)
#' 
#' Calculates the specified set of graphlet-based degree distributions from an 
#' indexed edge list using the ORCA fast graphlet orbit counting package.
#' @param indexed_edges A 2 x numEdges edgelist with vertices labelled with 
#' integer indices, with an optional "vertex_names" attribute
#' @param type Type of graphlet-based degree distributions: "orb4" counts all 
#' orbits for graphlets comprising up to 4 nodes; "orb5" counts all orbits for
#' graphlets comprising up to 5 nodes.
#' @return List of graphlet-based degree distributions, with each distribution
#' represented as a \code{dhist} discrete histogram object.
#' @export
gdd <- function(indexed_edges, type = "orb4") {
  orca_fn <- switch(type,
         orb4 = orca::count4,
         orb5 = orca::count5)
  orca_counts_to_graphlet_orbit_degree_distribution(orca_fn(indexed_edges))
}

#' Load all graphs in a directory and calculates their Graphlet-based Degree
#' Distributions (GDDs)
#' 
#' Loads graphs from all files matching the given pattern in the given directory,
#' converts them to indexed edge lists compatible with the ORCA fast orbit 
#' counting package and calculates the specified set of graphlet-based degree 
#' distributions usingthe ORCA package.
#' @param source_dir Path to graph directory
#' @param format Format of graph files
#' @param pattern Filename pattern to match graph files
#' @param type Type of graphlet-based degree distributions: "orb4" counts all 
#' orbits for graphlets comprising up to 4 nodes; "orb5" counts all orbits for
#' graphlets comprising up to 5 nodes.
#' @return A named list where each element contains a set of GDDs for a single 
#' graph from the source directory. Each set of GDDs is itself a named list,  
#' where each GDD element is a \code{dhist} discrete histogram object.
#' @export
gdd_for_all_graphs <- function(
  source_dir, format = "ncol", pattern = ".txt", type = "node4", 
  mc.cores = getOption("mc.cores", 2L)) {
  # Read graphs from source directory as ORCA-compatible edge lists
  edges <- read_all_graphs_as_orca_edge_lists(
    source_dir = source_dir, format = format, pattern = pattern)
  # Calculate specified GDDs for each graph
  parallel::mcmapply(gdd, edges, MoreArgs = list(type = type), 
                     SIMPLIFY = FALSE, mc.cores = mc.cores)
}
