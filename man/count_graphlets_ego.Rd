% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orca_interface.R
\name{count_graphlets_ego}
\alias{count_graphlets_ego}
\title{Ego-network graphlet counts}
\usage{
count_graphlets_ego(
  graph,
  max_graphlet_size = 4,
  neighbourhood_size = 2,
  min_ego_nodes = 3,
  min_ego_edges = 1,
  return_ego_networks = FALSE
)
}
\arguments{
\item{graph}{An undirected, simple graph as an \code{igraph} object.}

\item{max_graphlet_size}{Determines the maximum size of graphlets to count.
Only graphlets containing up to \code{max_graphlet_size} nodes will be
counted. Currently only size 4 (default) and 5 are supported.}

\item{neighbourhood_size}{The number of steps from the source node used to select the
neighboring nodes to be included in the source node ego-network. (Default 2).}

\item{min_ego_nodes}{Only ego networks with at least \code{min_ego_nodes}
nodes are returned. (Default 3).}

\item{min_ego_edges}{Only ego networks with at least \code{min_ego_edges}
edges are returned. (Default 1).}

\item{return_ego_networks}{If \code{TRUE}, return ego-networks alongside
graphlet counts to enable further processing.}
}
\value{
If \code{return_ego_networks = FALSE}, returns an RxC matrix
containing counts of each graphlet (columns, C) for each ego-network in the
input graph (rows, R). Columns are labelled with graphlet IDs and rows are
labelled with the ID of the central node in each ego-network (if nodes in the
input graph are labelled). If \code{return_ego_networks = TRUE}, returns a
list with the following elements:
\itemize{
  \item \code{graphlet_counts}: A matrix containing graphlet counts for each
  ego-network in the input graph as described above.
  \item \code{ego_networks}: The ego-networks of the query graph.
}
}
\description{
Calculates graphlet counts for the n-step ego-network of each node in a graph
}
