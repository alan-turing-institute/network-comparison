## ------------------------------------------------------------------------
# Load libraries
library("netdist")
library("purrr")

## ------------------------------------------------------------------------
# Set source directory for Virus PPI graph edge files
source_dir <- system.file(file.path("extdata", "VRPINS"), package = "netdist")

# Load query graphs
graph_1 <- read_simple_graph(file.path(source_dir, "EBV.txt"),
                             format = "ncol")

graph_2 <- read_simple_graph(file.path(source_dir, "ECL.txt"),
                             format = "ncol")


## ------------------------------------------------------------------------
# Maximum graphlet size to calculate counts and netdis statistic for.
max_graphlet_size <- 4

# Ego network neighbourhood size
neighbourhood_size <- 2

# Minimum size of ego networks to consider
min_ego_nodes <- 3
min_ego_edges <- 1

# Ego network density binning parameters
min_bin_count <- 5
num_bins <- 100

## ------------------------------------------------------------------------
# Get ego networks for query graphs and reference graph
ego_1 <- make_named_ego_graph(graph_1, 
                              order = neighbourhood_size, 
                              min_ego_nodes = min_ego_nodes, 
                              min_ego_edges = min_ego_edges)

ego_2 <- make_named_ego_graph(graph_2, 
                              order = neighbourhood_size, 
                              min_ego_nodes = min_ego_nodes, 
                              min_ego_edges = min_ego_edges)


## ------------------------------------------------------------------------
# Count graphlets for ego networks in query and reference graphs
graphlet_counts_1 <- ego_to_graphlet_counts(ego_1, max_graphlet_size = max_graphlet_size)
graphlet_counts_2 <- ego_to_graphlet_counts(ego_2, max_graphlet_size = max_graphlet_size)

## ------------------------------------------------------------------------
# rep(1, nrow(graphlet_counts)): list of ones as bin index, i.e. everything in same bin
mean_graphlet_counts_1 <- density_binned_counts(graphlet_counts_1,
                                                rep(1, nrow(graphlet_counts_1)))

mean_graphlet_counts_2 <- density_binned_counts(graphlet_counts_2,
                                                rep(1, nrow(graphlet_counts_2)))

bins <- c(0, 1)

## ------------------------------------------------------------------------
# Calculate expected graphlet counts for each ego network
exp_graphlet_counts_1 <- netdis_expected_graphlet_counts_per_ego(ego_1, 
                                                                 bins, 
                                                                 mean_graphlet_counts_1,
                                                                 max_graphlet_size,
                                                                 scale_fn = NULL)


exp_graphlet_counts_2 <- netdis_expected_graphlet_counts_per_ego(ego_2, 
                                                                 bins, 
                                                                 mean_graphlet_counts_2,
                                                                 max_graphlet_size,
                                                                 scale_fn = NULL)
# Centre graphlet counts by subtracting expected counts
centred_graphlet_counts_1 <- graphlet_counts_1 - exp_graphlet_counts_1

centred_graphlet_counts_2 <- graphlet_counts_2 - exp_graphlet_counts_2

## ------------------------------------------------------------------------
sum_graphlet_counts_1 <- colSums(centred_graphlet_counts_1)

sum_graphlet_counts_2 <- colSums(centred_graphlet_counts_2)

## ------------------------------------------------------------------------

netdis_result <- netdis_uptok(sum_graphlet_counts_1, 
                              sum_graphlet_counts_2, 
                              max_graphlet_size)

print(netdis_result)

