## ------------------------------------------------------------------------
# Load libraries
library("netdist")
library("purrr")

## ------------------------------------------------------------------------
# Maximum graphlet size to calculate counts for.
# We choose the specific graphlet size for the Netdis metric later.
max_graphlet_size = 4

## ------------------------------------------------------------------------
# Set source directory for Virus PPI graph edge files
source_dir <- system.file(file.path("extdata", "VRPINS"), package = "netdist")
# Load query graphs
graphs <- read_simple_graphs(source_dir, format = "ncol", pattern = "*")

## ------------------------------------------------------------------------

# Set ego network neighbourhood size
neighbourhood_size = 2
ego_networks <- purrr::map(graphs, make_named_ego_graph, 
                           order = neighbourhood_size)

## ------------------------------------------------------------------------
ego_graphlet_counts <- purrr::map_depth(ego_networks, 2, count_graphlets_for_graph,
                              max_graphlet_size = max_graphlet_size)

## ------------------------------------------------------------------------
# Load reference graph
file <- system.file(file.path("extdata", "random", "ER_1250_10_1"), 
                    package = "netdist")
ref_graph <- read_simple_graph(file, format = "ncol")
# Generate ego networks for reference graph
ref_ego_networks <- make_named_ego_graph(ref_graph, order = neighbourhood_size)
# Count graphlets for ego networks in reference graph
ref_ego_graphlet_counts <- purrr::map(ref_ego_networks, count_graphlets_for_graph,
                              max_graphlet_size = max_graphlet_size)

