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
#' connected, undirected and simple to ensure compatiblity with the ORCA
#' fast orbit calculation package:
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
#' Loads a graph from file as a connected, undirected, simple 
#' \code{igraph} object. 
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

#' Load all graphs in a directory, converting to ORCA compatible graphs
#' 
#' Loads graphs from all files matching the given pattern in the given directory
#' as connected, undirected, simple \code{igraph} objects. 
#' @param source_dir Path to graph directory
#' @param format Format of graph files
#' @param pattern Filename pattern to match graph files
#' @return A list of ORCA compatible graphs, with element names corresponding 
#' to the names of the files each graph was loaded from
#' @export
read_all_graphs_as_orca_graphs <- function (source_dir, format = "ncol", pattern = ".txt") {
  # Get list of all filenames in firectory that match the pattern
  file_names <- dir(source_dir, pattern = pattern)
  # Read graph data from each ".txt" file as an ORCA-compatible indexed edge list
  graphs <- purrr::map(file_names, function(file_name) {
    read_orca_graph(file = file.path(source_dir, file_name), format = "ncol")
  })
  # Name each graph with the name of the file it was read from
  attr(graphs, "names") <- file_names
  return(graphs)
}

#' ORCA vertex graphlet or orbit counts to graphlet-based histograms
#' 
#' Converts ORCA output (counts of each graphlet or orbit at each graph vertex) to 
#' a set of graphlet degree histograms (a histogram of counts across all graph 
#' vertices for each graphlet or orbit) 
#' @param orca_counts ORCA output: Counts of each graphlet or orbit 
#' (columns) at each graph vertex (rows)
#' @return Graphlet degree histograms: List of degree histograms for each 
#' graphlet or orbit
#' @export
orca_counts_to_graphlet_orbit_degree_distribution <- function(orca_counts) {
  apply(orca_counts, 2, dhist_from_obs)
}

#' Graphlet-based degree distributions (GDDs)
#' 
#' Generates graphlet-based degree distributions from \code{igraph} graph object,
#' using the ORCA fast graphlet orbit counting package.
#' @param graph A connected, undirected, simple graph as an \code{igraph} object. 
#' @param feature_type Type of graphlet-based feature to count: "graphlet"
#' counts the number of graphlets each node participates in; "orbit" calculates
#' the number of graphlet orbits each node participates in.
#' @param max_graphlet_size Determines the maximum size of graphlets to count. 
#' Only graphlets containing up to \code{max_graphlet_size} nodes will be counted.
#' @return List of graphlet-based degree distributions, with each distribution
#' represented as a \code{dhist} discrete histogram object.
#' @export
gdd <- function(graph, feature_type = 'orbit', max_graphlet_size = 4, 
                ego_neighbourhood_size = 0){
  if(ego_neighbourhood_size > 0) {
    if(feature_type != 'graphlet') {
      stop("Feature type not supported for ego-networks")
    } else {
      out <- count_graphlets_ego(graph, max_graphlet_size = max_graphlet_size, 
                                 neighbourhood_size = ego_neighbourhood_size)
    }
  } else  if(feature_type == "orbit") {
    out <- count_orbits(graph, max_graphlet_size = max_graphlet_size)
  } else if(feature_type == "graphlet") {
    out <- count_graphlets(graph, max_graphlet_size = max_graphlet_size)
  } else {
    stop('gdd: unrecognised feature_type')
  }
  orca_counts_to_graphlet_orbit_degree_distribution(out)
}

#' Calculate graphlet orbit counts
#' 
#' Calculates graphlet orbit counts from \code{igraph} graph object,
#' using the ORCA fast graphlet orbit counting package.
#' @param graph A connected, undirected, simple graph as an \code{igraph} object. 
#' @param max_graphlet_size Determines the maximum size of graphlets to count. 
#' Only graphlets containing up to \code{max_graphlet_size} nodes will be counted.
#' @return ORCA-format matrix containing counts of each graphlet
#' orbit (columns) at each vertex in the graph (rows).
#' @export
count_orbits <- function(graph, max_graphlet_size) {
    if(max_graphlet_size == 4) {
      orca_fn <- orca::count4
    } else if(max_graphlet_size == 5) {
      orca_fn <- orca::count5
    } else {
      stop("Unsupported maximum graphlet size")
    }
    indexed_edges <- graph_to_indexed_edges(graph)
    orbit_counts <- orca_fn(indexed_edges)
    rownames(orbit_counts) <- igraph::get.vertex.attribute(graph, name = "name")
    return(orbit_counts)
}

#' Calculate graphlet counts
#' 
#' Calculates graphlet counts from \code{igraph} graph object using the ORCA
#' fast graphlet orbit counting package by summing orbits over graphlets.
#' @param graph A connected, undirected, simple graph as an \code{igraph} object. 
#' @param max_graphlet_size Determines the maximum size of graphlets to count. 
#' Only graphlets containing up to \code{max_graphlet_size} nodes will be counted.
#' @return ORCA-format matrix containing counts of each graphlet (columns) at 
#' each vertex in the graph (rows).
#' @export
count_graphlets <- function(graph, max_graphlet_size) {
  orbit_counts <- count_orbits(graph, max_graphlet_size = max_graphlet_size)
  orbit_to_graphlet_counts(orbit_counts)
}

#' Ego-network graphlet counts
#' 
#' Calculates graphlet counts for the n-step ego-network of each node in a graph
#' @param graph A connected, undirected, simple graph as an \code{igraph} object. 
#' @param max_graphlet_size Determines the maximum size of graphlets to count. 
#' Only graphlets containing up to \code{max_graphlet_size} nodes will be counted.
#' @param neighbourhood_size The number of steps from the source node to include
#' nodes for each ego-network.
#' @param return_ego_networks If \code{TRUE}, return ego-networks alongside 
#' graphlet counts to enable further processing. 
#' @return If \code{return_ego_networks = FALSE}, returns an RxC matrix 
#' containing counts of each graphlet (columns, C) for each ego-network in the 
#' input graph (rows, R). Columns are labelled with graphlet IDs and rows are 
#' labelled with the ID of the central node in each ego-network (if nodes in the
#' input graph are labelled). If \code{return_ego_networks = TRUE}, returns a
#' list with the following elements:
#' \itemize{
#'   \item \code{graphlet_counts}: A matrix containing graphlet counts for each 
#'   ego-network in the input graph as described above.
#'   \item \code{ego_networks}: The ego-networks of the query graph.
#' }
#' @export
count_graphlets_ego <- function(graph, max_graphlet_size = 4, neighbourhood_size, 
                                return_ego_networks = FALSE) {
  # Extract ego network for each node in original graph, naming each ego network
  # in the list with the name of the node the ego network is generated for
  ego_networks <- make_named_ego_graph(graph, order = neighbourhood_size)
  # Generate graphlet counts for each node in each ego network (returns an ORCA
  # format graphlet count matrix for each ego network)
  ego_graphlet_counts_per_node <- purrr::map(ego_networks, count_graphlets, 
                                          max_graphlet_size = max_graphlet_size)
  # Sum graphlet counts across all nodes in each ego network
  ego_graphlet_counts <- purrr::map(ego_graphlet_counts_per_node, apply, MARGIN = 2, FUN = sum)
  # To ensure we only count each graphlet present in an ego network once, divide
  # the ego network graphlet counts by the number of nodes that contribute to 
  # each graphlet type
  ego_graphlet_counts <- purrr::map(ego_graphlet_counts, function(egc) {
    egc / graphlet_key(max_graphlet_size)$node_count})
  # Reshape the list of per node single row graphlet count matrices to a single
  # ORCA format graphlet count matrix with one row per node
  ego_graphlet_counts <- t(simplify2array(ego_graphlet_counts))
  # Return either graphlet counts, or graphlet counts and ego_networks
  if(return_ego_networks) {
    return(list(graphlet_counts = ego_graphlet_counts, ego_networks = ego_networks))
  } else {
    return(ego_graphlet_counts)
  }
}

#' Get ego-networks for a graph as a named list
#'
#' Simple wrapper for the \code{igraph::make_ego_graph} function that names
#' each ego-network in the returned list with the name of the node in the
#' original graph that the ego-network was generated from
#' @export
make_named_ego_graph <- function(graph, order, ...) {
  ego_networks <- igraph::make_ego_graph(graph, order, ...)
  names(ego_networks) <- igraph::V(graph)$name
  ego_networks
}

#' Orbit to graphlet counts
#' 
#' Converts graphlet orbit counts at each vertex to graphlet counts at each 
#' vertex by summing over all orbits contained within each graphlet
#' @param orbit_counts ORCA-format matrix containing counts of each graphlet
#' orbit (columns) at each vertex in the graph (rows)
#' @return An ORCA-style matrix containing counts of each graphlet (columns) at
#' each vertex in the graph (rows)
#' @export
orbit_to_graphlet_counts <- function(orbit_counts) {
  num_orbits <- dim(orbit_counts)[2]
  # Indexes to select the orbit(s) that comprise each graphlet. Note that we 
  # define these in the zero-based indexing used in journal papers, but 
  # need to add 1 to convert to the 1-based indexing used by R
  if(num_orbits == 15) {
    # Orbits for graphlets comprising up to 4 nodes
    max_nodes <- 4
    orbit_to_graphlet_map <- 
      purrr::map(list(0, 1:2, 3, 4:5, 6:7, 8, 9:11, 12:13, 14), 
                 function(indexes){ indexes + 1})
  } else if(num_orbits == 73) {
    # Orbits for graphlets comprising up to 5 nodes
    max_nodes <- 5
    orbit_to_graphlet_map <- 
      purrr::map(list(0, 1:2, 3, 4:5, 6:7, 8, 9:11, 12:13, 14, 15:17, 18:21, 
                      22:23, 24:26, 27:30, 31:33, 34, 35:38, 39:42, 43:44, 
                      45:48, 49:50, 51:53, 54:55, 56:58, 59:61, 62:64, 
                      65:67, 68:69, 70:71, 72), 
                 function(indexes){ indexes + 1})
  } else {
    stop(("Unsupported number of orbits"))
  }
  # Sum counts across orbits in graphlets
  graphlet_counts <- sapply(orbit_to_graphlet_map, function(indexes){
    rowSums(orbit_counts[,indexes, drop = FALSE])})
  # Add graphlet names
  colnames(graphlet_counts) <- graphlet_key(max_nodes)$id
  return(graphlet_counts)
}

#' Graphlet key
#' 
#' Metdata about graphlet groups.
#' @param max_graphlet_size Maximum number of nodes graphlets can contain
#' @return Metadata list with the following named fields:
#' \itemize{
#'   \item \code{max_nodes}: Maximum number of nodes graphlets can contain
#'   \item \code{id}: ID of each graphlet in format Gn, where n is in range 0 to 
#'  num_graphlets
#'   \item \code{node_count}: Number of nodes contained within each graphlet
#' }
#' @export
graphlet_key <- function(max_graphlet_size) {
  if(max_graphlet_size == 4) {
    node_count <- c(2, rep(3,2), rep(4,6))
  } else if (max_graphlet_size == 5) {
    node_count <- c(2, rep(3,2), rep(4,6), rep(5, 21))
  } else {
    stop("Unsupported maximum graphlet size")
  }
  max_node_index <- length(node_count)-1
  id <- purrr::simplify(purrr::map(0:max_node_index, function(index) {
    paste('G', index, sep = "")}))
  name <- 
  return(list(max_nodes = max_graphlet_size, id = id, node_count = node_count))
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
  source_dir, format = "ncol", pattern = ".txt", feature_type = "orbit", 
  max_graphlet_size = 4, mc.cores = getOption("mc.cores", 2L)) {
  # Read graphs from source directory as ORCA-compatible edge lists
  graphs <- read_all_graphs_as_orca_graphs(
    source_dir = source_dir, format = format, pattern = pattern)
  # Calculate specified GDDs for each graph
  parallel::mcmapply(gdd, graphs, MoreArgs = 
                       list(feature_type = feature_type, 
                            max_graphlet_size = max_graphlet_size), 
                     SIMPLIFY = FALSE, mc.cores = mc.cores)
}
