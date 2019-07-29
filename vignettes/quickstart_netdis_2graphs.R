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

# Load reference graph
# JACK - need to deal with case where ref graph not used.
ref_path <- system.file(file.path("extdata", "random", "ER_1250_10_1"), 
                        package = "netdist")
ref_graph <- read_simple_graph(ref_path, format = "ncol")

## ------------------------------------------------------------------------
# Maximum graphlet size to calculate counts and netdis statistic for.
max_graphlet_size = 4

# Ego network neighbourhood size
neighbourhood_size = 2

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

ego_ref <- make_named_ego_graph(ref_graph, 
                                order = neighbourhood_size, 
                                min_ego_nodes = min_ego_nodes, 
                                min_ego_edges = min_ego_edges)

## ------------------------------------------------------------------------
# Count graphlets for ego networks in query and reference graphs
graphlet_counts_1 <- ego_to_graphlet_counts(ego_1, max_graphlet_size = max_graphlet_size)
graphlet_counts_2 <- ego_to_graphlet_counts(ego_2, max_graphlet_size = max_graphlet_size)

graphlet_counts_ref <- ego_to_graphlet_counts(ego_ref, max_graphlet_size = max_graphlet_size)

## ------------------------------------------------------------------------

# Scale ego-network graphlet counts by dividing by total number of k-tuples in
# ego-network (where k is graphlet size)
scaled_graphlet_counts_ref <- scale_graphlet_counts_ego(ego_ref, 
                                                        graphlet_counts_ref, 
                                                        max_graphlet_size)

# Get ego-network densities
densities_ref <- ego_network_density(ego_ref)

# Adaptively bin ref ego-network densities
binned_densities <- binned_densities_adaptive(densities_ref, 
                                              min_counts_per_interval = min_bin_count, 
                                              num_intervals = num_bins)

ref_ego_density_bins <- binned_densities$breaks

# Average ref graphlet counts across density bins
ref_binned_graphlet_counts <- mean_density_binned_graphlet_counts(
                                  scaled_graphlet_counts_ref, 
                                  binned_densities$interval_indexes)
  

## ------------------------------------------------------------------------
# Calculate expected graphlet counts (using ref graph ego network density bins)
exp_graphlet_counts_1 <- netdis_expected_graphlet_counts_per_ego(ego_1, 
                                                                 max_graphlet_size,
                                                                 ref_ego_density_bins, 
                                                                 ref_binned_graphlet_counts)


exp_graphlet_counts_2 <- netdis_expected_graphlet_counts_per_ego(ego_2, 
                                                                 max_graphlet_size,
                                                                 ref_ego_density_bins, 
                                                                 ref_binned_graphlet_counts)

# Centre graphlet counts by subtracting expected counts
centred_graphlet_counts_1 <- graphlet_counts_1 - exp_graphlet_counts_1

centred_graphlet_counts_2 <- graphlet_counts_2 - exp_graphlet_counts_2

## ------------------------------------------------------------------------
sum_graphlet_counts_1 <- colSums(centred_graphlet_counts_1)

sum_graphlet_counts_2 <- colSums(centred_graphlet_counts_2)

## ------------------------------------------------------------------------
netdis_uptok(sum_graphlet_counts_1, 
             sum_graphlet_counts_2, 
             max_graphlet_size)


## ------------------------------------------------------------------------


