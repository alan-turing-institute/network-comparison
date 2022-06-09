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
  graph <- igraph::set.vertex.attribute(
    graph,
    name = "name",
    value = attr(indexed_edges, "vertex_names")
  )
  return(graph)
}

#' Read all graphs in a directory, simplifying as requested
#'
#' Reads graph data from all files in a directory matching the specified
#' filename pattern. From each file, an a igraph graph object is constructed
#' and the requested subset of the following simplifications is made in the
#' following order:
#'   1. Makes the graph undirected
#'   2. Removes loops (where both endpoints of an edge are the same vertex)
#'   3. Removes multiple edges (i.e. ensuring only one edge exists for each
#'      pair of endpoints)
#'   4. Removes isolated vertices (i.e. vertices with no edges after the
#'      previous alterations)
#' @param source_dir Path to directory containing files with graph data
#' @param format Format of graph data. Any format supported by
#' \code{igraph::read_graph} can be used.
#' @param pattern Pattern to use to filter filenames. Any pattern supported by
#' \code{dir} can be used.
#' @param as_undirected If TRUE make graph edges undirected
#' @param remove_loops If TRUE, remove edgeds that connect a vertex to itself
#' @param remove_multiple If TRUE remove multiple edges connencting the same
#' pair of vertices
#' @param remove_isolates If TRUE, remove vertices with no edges after the
#' previous alterations have been made
#' @return A named list of simplified igraph graph object, with the name of each
#' graph set to the name of the file it was read from.
#' @examples
#' # Set source directory for Virus protein-protein interaction edge files
#' # stored in the netdist package.
#' source_dir <- system.file(
#'   file.path("extdata", "VRPINS"),
#'   package = "netdist"
#' )
#' print(source_dir)
#' # Load query graphs as igraph objects
#' graph_1 <- read_simple_graph(
#'   file.path(source_dir, "EBV.txt"),
#'   format = "ncol"
#' )
#' graph_1
#' @export
read_simple_graphs <- function(source_dir,
                               format = "ncol",
                               pattern = "*",
                               as_undirected = TRUE,
                               remove_loops = TRUE,
                               remove_multiple = TRUE,
                               remove_isolates = TRUE) {
  # Get list of all filenames in directory that match the pattern
  file_names <- dir(source_dir, pattern = pattern)
  # Read graph data from each matched file as an igraph format graph,
  # simplifying as requested
  graphs <- purrr::map(
    file_names,
    function(file_name) {
      read_simple_graph(
        file = file.path(source_dir, file_name),
        format = format,
        as_undirected = as_undirected,
        remove_loops = remove_loops,
        remove_multiple = remove_multiple,
        remove_isolates = remove_isolates
      )
    }
  )

  # Name each graph with the name of the file it was read from (with any
  # extension moved)
  names <- purrr::simplify(
    purrr::map(
      strsplit(file_names, "\\."),
      function(s) {
        if (length(s) == 1) {
          s
        } else {
          paste(utils::head(s, -1), collapse = ".")
        }
      }
    )
  )
  attr(graphs, "names") <- names
  return(graphs)
}

#' Read a graph from file, simplifying as requested
#'
#' Reads graph data from file, constructing an a igraph graph object, making the
#' requested subset of the following simplifications in the following order:
#'   1. Makes the graph undirected
#'   2. Removes loops (where both endpoints of an edge are the same vertex)
#'   3. Removes multiple edges (i.e. ensuring only one edge exists for each
#'      pair of endpoints)
#'   4. Removes isolated vertices (i.e. vertices with no edges after the
#'      previous alterations).
#' @param file Path to file containing graph data
#' @param format Format of graph data. All formats supported by
#' \code{igraph::read_graph} are supported.
#' @param as_undirected If TRUE make graph edges undirected
#' @param remove_loops If TRUE, remove edgeds that connect a vertex to itself
#' @param remove_multiple If TRUE remove multiple edges connencting the same
#' pair of vertices
#' @param remove_isolates If TRUE, remove vertices with no edges after the
#' previous alterations have been made
#' @return A simplified igraph graph object
#' @export
read_simple_graph <- function(file, format, as_undirected = TRUE,
                              remove_loops = TRUE, remove_multiple = TRUE,
                              remove_isolates = TRUE) {
  # Read graph from file. NOTE: igraph only supported the "directed" argument
  # for some formats, but passes it to formats that don't support it, which
  # then throw an error
  if (format %in% c("edgelist", "ncol", "lgl", "dimacs", "dl")) {
    graph <- igraph::read_graph(file = file, format = format, directed = TRUE)
  } else {
    graph <- igraph::read_graph(file = file, format = format)
  }
  # Perform any requested simplifications
  simplify_graph(graph,
    as_undirected = as_undirected,
    remove_loops = remove_loops, remove_multiple = remove_multiple,
    remove_isolates = remove_isolates
  )
}

#' Simplify an igraph
#'
#' Takes a igraph graph object and makes the requested subset of the following
#' simplifications in the following order:
#'   1. Makes the graph undirected
#'   2. Removes loops (where both endpoints of an edge are the same vertex)
#'   3. Removes multiple edges (i.e. ensuring only one edge exists for each
#'      pair of endpoints)
#'   4. Removes isolated vertices (i.e. vertices with no edges after the
#'      previous alterations)
#' @param graph An graph or list of graphs in igraph format
#' @param as_undirected If TRUE make graph edges undirected
#' @param remove_loops If TRUE, remove edgeds that connect a vertex to itself
#' @param remove_multiple If TRUE remove multiple edges connencting the same
#' pair of vertices
#' @param remove_isolates If TRUE, remove vertices with no edges after the
#' previous alterations have been made
#' @return A simplified igraph graph object
#' @export
simplify_graph <- function(graph, as_undirected = TRUE, remove_loops = TRUE,
                           remove_multiple = TRUE, remove_isolates = TRUE) {
  if (as_undirected) {
    # Ensure graph is undirected
    graph <- igraph::as.undirected(graph, mode = "each")
  }
  if (remove_loops || remove_multiple) {
    # Remove loops (where both endpoints of an edge are the same vertex) and
    # multiple edges (where two edges have the same endpoints [in the same order
    # for directed graphs])
    graph <- igraph::simplify(graph,
      remove.loops = remove_loops,
      remove.multiple = remove_multiple
    )
  }
  if (remove_isolates) {
    # Remove vertices that have no edges connecting them to other vertices
    # NOTE: Vertices that only connect to themselves will only be removed if
    # their self-connecting edges have been removed by setting remove_loops to
    # TRUE
    isolated_vertex_indices <- (igraph::degree(graph) == 0)
    graph <- igraph::delete.vertices(graph, isolated_vertex_indices)
  }
  return(graph)
}

#' Convert a matrix of node level features to a "discrete histogram" for
#' each feature.
#'
#' Converts a matrix of node level features (e.g. for example counts
#' of multiple graphlets or orbits at each node) to
#' a set of histogram like objects (observed frequency distribution of each
#' feature/column)
#' @param features_matrix A matrix whose rows represent nodes and whose columns
#' represent different node level features. This means that entry ij provides
#' the value of feature j for node i.
#' @return Feature histograms: List of "discrete histograms" for each feature
#' @export
graph_features_to_histograms <- function(features_matrix) {
  apply(features_matrix, 2, dhist_from_obs)
}

#' Graphlet-based degree distributions (GDDs)
#'
#' Short-cut function to create graphlet-based degree distributions from
#' \code{igraph} graph object using the ORCA fast graphlet orbit counting
#' package.
#' @param graph A connected, undirected, simple graph as an \code{igraph} object
#' @param feature_type Type of graphlet-based feature to count: "graphlet"
#' counts the number of graphlets each node participates in; "orbit" calculates
#' the number of graphlet orbits each node participates in.
#' @param max_graphlet_size Determines the maximum size of graphlets to count.
#' Only graphlets containing up to \code{max_graphlet_size} nodes will be
#' counted. Currently only size 4 and 5 are supported.
#' @param ego_neighbourhood_size The number of steps from the source node used
#' to select the
#' neighboring nodes to be included in the source node ego-network.
#' @return List of graphlet-based degree distributions, with each distribution
#' represented as a \code{dhist} discrete histogram object.
#' @export
gdd <- function(graph, feature_type = "orbit", max_graphlet_size = 4,
                ego_neighbourhood_size = 0) {
  graph <- simplify_graph(graph)
  if (ego_neighbourhood_size > 0) {
    if (feature_type != "graphlet") {
      stop("Feature type not supported for ego-networks")
    } else {
      out <- count_graphlets_ego(graph,
        max_graphlet_size = max_graphlet_size,
        neighbourhood_size = ego_neighbourhood_size
      )
    }
  } else if (feature_type == "orbit") {
    out <- count_orbits_per_node(graph, max_graphlet_size = max_graphlet_size)
  } else if (feature_type == "graphlet") {
    out <- count_graphlets_per_node(graph,
      max_graphlet_size = max_graphlet_size
    )
  } else {
    stop("gdd: unrecognised feature_type")
  }
  graph_features_to_histograms(out)
}

#' Count graphlet orbits for each node in a graph
#'
#' Calculates graphlet orbit counts for each node in an \code{igraph} graph
#' object, using the \code{orca} fast graphlet orbit counting package.
#' @param graph A undirected, simple graph as an \code{igraph} object.
#' @param max_graphlet_size Determines the maximum size of graphlets to count.
#' Only graphlets containing up to \code{max_graphlet_size} nodes will be
#' counted. Currently only size 4 and 5 are supported.
#' @return ORCA-format matrix containing counts of each graphlet
#' orbit (columns) at each node in the graph (rows).
#' @export
count_orbits_per_node <- function(graph, max_graphlet_size) {
  if (max_graphlet_size == 4) {
    orca_fn <- orca::count4
  } else if (max_graphlet_size == 5) {
    orca_fn <- orca::count5
  } else {
    stop("Unsupported maximum graphlet size")
  }
  indexed_edges <- graph_to_indexed_edges(graph)
  num_edges <- dim(indexed_edges)[[1]]
  if (num_edges >= 1) {
    orbit_counts <- orca_fn(indexed_edges)
  } else {
    # ORCA functions expect at least one edge, so handle this case separately
    # and manually construct empty orbit count matrix
    orbit_ids <- orbit_key(max_graphlet_size)$id
    num_orbits <- length(orbit_ids)
    num_nodes <- igraph::vcount(graph)
    orbit_counts <- matrix(0, nrow = num_nodes, ncol = num_orbits)
    colnames(orbit_counts) <- orbit_ids
  }
  rownames(orbit_counts) <- igraph::get.vertex.attribute(graph, name = "name")
  return(orbit_counts)
}

#' Count graphlets for each node in a graph
#'
#' Calculates graphlet counts for each node in an \code{igraph} graph object,
#' using the ORCA fast graphlet orbit counting package. by summing orbits over
#' graphlets.
#' @param graph A connected, undirected, simple graph as an \code{igraph} object
#' @param max_graphlet_size Determines the maximum size of graphlets to count.
#' Only graphlets containing up to \code{max_graphlet_size} nodes will be
#' counted. Currently only size 4 and 5 are supported.
#' @return ORCA-format matrix containing counts of each graphlet (columns) at
#' each node in the graph (rows).
#' @export
count_graphlets_per_node <- function(graph, max_graphlet_size) {
  orbit_counts <- count_orbits_per_node(graph,
    max_graphlet_size = max_graphlet_size
  )
  orbit_to_graphlet_counts(orbit_counts)
}

#' Count total number of graphlets in a graph
#'
#' Calculates total graphlet counts for a \code{igraph} graph object using the
#' ORCA fast graphlet orbit counting package. Per-node graphlet counts are
#' calculated by summing orbits over graphlets. These are then divided by the
#' number of nodes comprising each graphlet to avoid counting the same graphlet
#' multiple times.
#' @param graph A connected, undirected, simple graph as an \code{igraph} object
#' @param max_graphlet_size Determines the maximum size of graphlets to count.
#' Only graphlets containing up to \code{max_graphlet_size} nodes will be
#' counted. Currently only size 4 and 5 are supported.
#' @return Vector containing counts of each graphlet for the graph.
#' @export
count_graphlets_for_graph <- function(graph, max_graphlet_size) {
  node_counts <- count_graphlets_per_node(graph, max_graphlet_size)
  # Sum graphlet counts over all nodes (rows)
  total_counts <- colSums(node_counts)
  # To ensure we only count each graphlet present in an ego network once, divide
  # the graphlet counts by the number of nodes that contribute to
  # each graphlet type
  nodes_per_graphlet <- graphlet_key(max_graphlet_size)$node_count
  total_counts <- total_counts / nodes_per_graphlet

  # add overall graph node count to total_counts
  N <- igraph::vcount(graph) # nolint: object_name_linter.
  total_counts <- c(N = N, total_counts)
  total_counts
}

#' Ego-network graphlet counts
#'
#' Calculates graphlet counts for the n-step ego-network of each node in a graph
#' @param graph An undirected, simple graph as an \code{igraph} object.
#' @param max_graphlet_size Determines the maximum size of graphlets to count.
#' Only graphlets containing up to \code{max_graphlet_size} nodes will be
#' counted. Currently only size 4 (default) and 5 are supported.
#' @param neighbourhood_size The number of steps from the source node used to
#' select the neighboring nodes to be included in the source node ego-network.
#' (Default 2).
#' @param min_ego_nodes Only ego networks with at least \code{min_ego_nodes}
#' nodes are returned. (Default 3).
#' @param min_ego_edges Only ego networks with at least \code{min_ego_edges}
#' edges are returned. (Default 1).
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
count_graphlets_ego <- function(graph,
                                max_graphlet_size = 4,
                                neighbourhood_size = 2,
                                min_ego_nodes = 3,
                                min_ego_edges = 1,
                                return_ego_networks = FALSE) {
  # Extract ego network for each node in original graph, naming each ego network
  # in the list with the name of the node the ego network is generated for
  ego_networks <- make_named_ego_graph(graph,
    order = neighbourhood_size,
    min_ego_nodes = min_ego_nodes,
    min_ego_edges = min_ego_edges
  )

  # Generate graphlet counts for each node in each ego network
  ego_graphlet_counts <- ego_to_graphlet_counts(ego_networks, max_graphlet_size)

  # Return either graphlet counts, or graphlet counts and ego_networks
  if (return_ego_networks) {
    return(list(
      graphlet_counts = ego_graphlet_counts,
      ego_networks = ego_networks
    ))
  } else {
    return(ego_graphlet_counts)
  }
}

#' ego_to_graphlet_counts
#'
#' Calculates graphlet counts for previously generated ego networks.
#' @param ego_networks Named list of ego networks for a graph.
#' @param max_graphlet_size Determines the maximum size of graphlets to count.
#' Only graphlets containing up to \code{max_graphlet_size} nodes will be
#' counted. Currently only size 4 and 5 are supported.
#' @return returns an RxC matrix
#' containing counts of each graphlet (columns, C) for each ego-network
#' (rows, R).
#' Columns are labelled with graphlet IDs and rows are
#' labelled with the ID of the central node in each ego-network.
#' @export
ego_to_graphlet_counts <- function(ego_networks, max_graphlet_size = 4) {
  # Generate graphlet counts for each node in each ego network (returns an ORCA
  # format graphlet count matrix for each ego network)
  ego_graphlet_counts <- purrr::map(ego_networks, count_graphlets_for_graph,
    max_graphlet_size = max_graphlet_size
  )

  # Reshape the list of per node single row graphlet count matrices to a single
  # ORCA format graphlet count matrix with one row per node
  ego_graphlet_counts <- t(simplify2array(ego_graphlet_counts))

  # Return graphlet counts
  return(ego_graphlet_counts)
}

#' Get ego-networks for a graph as a named list
#'
#' Simple wrapper for the \code{igraph::make_ego_graph} function that names
#' each ego-network in the returned list with the name of the node in the
#' original graph that the ego-network was generated from
#' @param graph An \code{igraph} object
#' @param order The number of steps from the source node to include
#' nodes for each ego-network.
#' @param min_ego_nodes Only ego networks with at least \code{min_ego_nodes}
#' nodes are returned.
#' @param min_ego_edges Only ego networks with at least \code{min_ego_edges}
#' edges are returned.
#' @param ... Additional parameters to be passed to the underlying
#' \code{igraph::make_ego_graph} function used.
#' @export
make_named_ego_graph <- function(graph, order, min_ego_nodes = 3,
                                 min_ego_edges = 1, ...) {
  ego_networks <- igraph::make_ego_graph(graph, order, ...)
  names(ego_networks) <- igraph::V(graph)$name

  # Drop ego-networks that don't have the minimum number of nodes or edges
  drop_index <- purrr::simplify(purrr::map(ego_networks, function(g) {
    (igraph::vcount(g) < min_ego_nodes) | (igraph::ecount(g) < min_ego_edges)
  }))
  ego_networks <- ego_networks[!drop_index]

  return(ego_networks)
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
  if (num_orbits == 15) {
    # Orbits for graphlets comprising up to 4 nodes
    max_nodes <- 4
    orbit_to_graphlet_map <-
      purrr::map(
        list(0, 1:2, 3, 4:5, 6:7, 8, 9:11, 12:13, 14),
        function(indexes) {
          indexes + 1
        }
      )
  } else if (num_orbits == 73) {
    # Orbits for graphlets comprising up to 5 nodes
    max_nodes <- 5
    orbit_to_graphlet_map <-
      purrr::map(
        list(
          0, 1:2, 3, 4:5, 6:7, 8, 9:11, 12:13, 14, 15:17, 18:21,
          22:23, 24:26, 27:30, 31:33, 34, 35:38, 39:42, 43:44,
          45:48, 49:50, 51:53, 54:55, 56:58, 59:61, 62:64,
          65:67, 68:69, 70:71, 72
        ),
        function(indexes) {
          indexes + 1
        }
      )
  } else {
    stop(("Unsupported number of orbits"))
  }
  # Sum counts across orbits in graphlets
  graphlet_counts <- sapply(orbit_to_graphlet_map, function(indexes) {
    rowSums(orbit_counts[, indexes, drop = FALSE])
  })
  if (dim(orbit_counts)[[1]] == 1) {
    # If orbit counts has only a single row, sapply returns a vector
    # rather than a matrix, so convert to a matrix by adding dim
    dim(graphlet_counts) <- c(1, length(graphlet_counts))
  }
  # Add graphlet names
  colnames(graphlet_counts) <- graphlet_key(max_nodes)$id
  return(graphlet_counts)
}

#' Graphlet key
#'
#' Metdata about graphlet groups.
#' @param max_graphlet_size Maximum number of nodes graphlets can contain.
#' Currently only size 2 to 5 are supported.
#' @return Metadata list with the following named fields:
#' \itemize{
#'   \item \code{max_nodes}: Maximum number of nodes graphlets can contain
#'   \item \code{id}: ID of each graphlet in format Gn, where n is in range 0 to
#'  num_graphlets
#'   \item \code{node_count}: Number of nodes contained within each graphlet
#' }
#' @export
graphlet_key <- function(max_graphlet_size) {
  if (max_graphlet_size == 2) {
    node_count <- c(2)
  } else if (max_graphlet_size == 3) {
    node_count <- c(2, rep(3, 2))
  } else if (max_graphlet_size == 4) {
    node_count <- c(2, rep(3, 2), rep(4, 6))
  } else if (max_graphlet_size == 5) {
    node_count <- c(2, rep(3, 2), rep(4, 6), rep(5, 21))
  } else {
    stop("Unsupported maximum graphlet size")
  }
  max_node_index <- length(node_count) - 1
  id <- purrr::simplify(purrr::map(0:max_node_index, function(index) {
    paste("G", index, sep = "")
  }))
  name <-
    return(list(
      max_nodes = max_graphlet_size,
      id = id,
      node_count = node_count
    ))
}

#' Orbit key
#'
#' Metdata about orbit groups.
#' @param max_graphlet_size Maximum number of nodes graphlets can contain.
#' Currently only size 2 to 5 are supported.
#' @return Metadata list with the following named fields:
#' \itemize{
#'   \item \code{max_nodes}: Maximum number of nodes graphlets can contain
#'   \item \code{id}: ID of each graphlet in format On, where n is in range 0 to
#'  num_orbits
#'   \item \code{node_count}: Number of nodes contained within each graphlet
#' }
#' @export
orbit_key <- function(max_graphlet_size) {
  if (max_graphlet_size == 2) {
    node_count <- c(2)
  } else if (max_graphlet_size == 3) {
    node_count <- c(2, rep(3, 3))
  } else if (max_graphlet_size == 4) {
    node_count <- c(2, rep(3, 3), rep(4, 11))
  } else if (max_graphlet_size == 5) {
    node_count <- c(2, rep(3, 3), rep(4, 11), rep(5, 58))
  } else {
    stop("Unsupported maximum graphlet size")
  }
  max_node_index <- length(node_count) - 1
  id <- purrr::simplify(purrr::map(0:max_node_index, function(index) {
    paste("O", index, sep = "")
  }))
  name <-
    return(list(
      max_nodes = max_graphlet_size,
      id = id,
      node_count = node_count
    ))
}

#' Graphlet IDs for size
#'
#' List IDs for all graphlets of a specified size
#' @param graphlet_size Graphlet size (i.e. number of nodes in the graphlet)
#' @return A vector containing the IDs of all graphlets made up of the specified
#' number of nodes
#' @export
graphlet_ids_for_size <- function(graphlet_size) {
  graphlet_key <- graphlet_key(graphlet_size)
  graphlet_key$id[graphlet_key$node_count == graphlet_size]
}

#' Load all graphs in a directory and calculates their Graphlet-based Degree
#' Distributions (GDDs)
#'
#' Loads graphs from all files matching the given pattern in the given
#' directory, converts them to indexed edge lists compatible with the ORCA fast
#' orbit counting package and calculates the specified set of graphlet-based
#' degree distributions usingthe ORCA package.
#' @param source_dir Path to graph directory
#' @param format Format of graph files
#' @param pattern Filename pattern to match graph files
#' @param feature_type Type of graphlet-based degree distributions. Can be
#' \code{graphlet} to count graphlets or \code{orbit} to count orbits.
#' @return A named list where each element contains a set of GDDs for a single
#' @param max_graphlet_size Maximum size of graphlets to use when generating
#' GDD. Currently only size 4 and 5 are supported.
#' @param ego_neighbourhood_size The number of steps from the source node used
#' to select the neighboring nodes to be included in the source node
#' ego-network. If set to 0, ego-networks will not be used.
#' @param  mc.cores Number of cores to use for parallel processing. Defaults to
#' the \code{mc.cores} option set in the R environment.
#' @return A named list where each element contains a set of GDDs for a single
#' graph from the source directory. Each set of GDDs is itself a named list,
#' where each GDD element is a \code{dhist} discrete histogram object.
#' @export
gdd_for_all_graphs <- function(source_dir,
                               format = "ncol",
                               pattern = ".txt",
                               feature_type = "orbit",
                               max_graphlet_size = 4,
                               ego_neighbourhood_size = 0,
                               mc.cores = getOption("mc.cores", 2L)) { # nolint: object_name_linter.
  # Create function to read graph from file and generate GDD
  graphs <- read_simple_graphs(
    source_dir = source_dir, format = format, pattern = pattern
  )
  # Calculate specified GDDs for each graph
  # NOTE: mcapply only works on unix-like systems with system level forking
  # capability. This means it will work on Linux and OSX, but not Windows.
  # For now, we just revert to single threaded operation on Windows
  # TODO: Look into using the parLappy function on Windows
  if (.Platform$OS.type != "unix") {
    # Force cores to 1 if system is not unix-like as it will not support
    # forking
    mc.cores <- 1 # nolint: object_name_linter.
  }
  parallel::mcmapply(gdd, graphs,
    MoreArgs =
      list(
        feature_type = feature_type,
        max_graphlet_size = max_graphlet_size,
        ego_neighbourhood_size = ego_neighbourhood_size
      ),
    SIMPLIFY = FALSE, mc.cores = mc.cores
  )
}

#' Generate a cross-comparison specification
#'
#' Creates a cross-comparison matrix with pair-wise combinations
#' of elements from the provided list.
#' @param named_list A named list of items for which an exhaustive pair-wise
#' cross-comparison is required.
#' @param how How to generate pair-wise combinations. Either "many-to-many"
#' (default) which generates all possible pair-wise combinations, or
#' "one-to-many" which generates all combinations between the first element
#' in named_list and the rest of the elements only.
#' @return A matrix with one row for each possible pair-wise combination
#' of elements from the provided named list. The first and second columns
#' contain the names of the elements in the pair and the third and fourth
#' columns contain the indexes of these elements in the provided list.
#' @export
cross_comparison_spec <- function(named_list, how = "many-to-many") {
  if (how == "one-to-many") {
    indexes <- data.frame(
      rep(1, length(named_list) - 1),
      2:length(named_list)
    )
  } else {
    indexes <- as.data.frame(t(utils::combn(1:length(named_list), 2)))
  }

  names <- as.data.frame(cbind(
    names(named_list)[indexes[, 1]],
    names(named_list)[indexes[, 2]]
  ))
  spec <- cbind(names, indexes)
  colnames(spec) <- c("name_a", "name_b", "index_a", "index_b")
  return(spec)
}

#' Convert a pair-wise cross-comparison into a matrix format
#'
#' Converts a pair-wise cross-comparison into a matrix format
#' @param measure A list of pair-wise comparison measiures
#' @param cross_comparison_spec A cross-comparison specification generated
#' using \code{cross_comparison_spec}
#' @return A square symmetric matrix with a zero diagonal, with elements
#' Cij and Cji populated from the element from \code{measure} corresponding to
#' the row of \code{cross_comparison_spec} with \code{index_a = i} and
#' \code{index_b = j}
#' @export
cross_comp_to_matrix <- function(measure, cross_comparison_spec) {
  num_items <- max(c(
    cross_comparison_spec$index_a,
    cross_comparison_spec$index_b
  ))
  out <- matrix(data = 0, nrow = num_items, ncol = num_items)
  out[cbind(
    cross_comparison_spec$index_a,
    cross_comparison_spec$index_b
  )] <- measure
  out[cbind(
    cross_comparison_spec$index_b,
    cross_comparison_spec$index_a
  )] <- measure
  row_labels <- rep("<MISSING>", num_items)
  row_labels[cross_comparison_spec$index_a] <- as.character(
    cross_comparison_spec$name_a
  )
  row_labels[cross_comparison_spec$index_b] <- as.character(
    cross_comparison_spec$name_b
  )
  rownames(out) <- row_labels
  col_labels <- rep("<MISSING>", num_items)
  col_labels[cross_comparison_spec$index_a] <- as.character(
    cross_comparison_spec$name_a
  )
  col_labels[cross_comparison_spec$index_b] <- as.character(
    cross_comparison_spec$name_b
  )
  colnames(out) <- col_labels
  return(out)
}
