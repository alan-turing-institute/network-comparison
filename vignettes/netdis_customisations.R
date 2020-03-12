## ------------------------------------------------------------------------
# Load libraries
library("netdist")
library("purrr")

## ------------------------------------------------------------------------
# Maximum graphlet size to calculate counts and netdis statistic for.
max_graphlet_size <- 4

# Ego network neighbourhood size
neighbourhood_size <- 2

# Minimum size of ego networks to consider
min_ego_nodes <- 3
min_ego_edges <- 1

# Reference graph
ref_path <- system.file(file.path("extdata", "random", "ER_1250_10_1"), 
                        package = "netdist")
ref_graph <- read_simple_graph(ref_path, format = "ncol")


## ------------------------------------------------------------------------
source_dir <- system.file(file.path("extdata", "VRPINS"), package = "netdist")
graphs <- read_simple_graphs(source_dir, format = "ncol", pattern = "*")

## ------------------------------------------------------------------------

# Calculate netdis statistics
results <- netdis_many_to_many(graphs,
                               ref_graph,
                               max_graphlet_size = max_graphlet_size,
                               neighbourhood_size = neighbourhood_size,
                               min_ego_nodes = min_ego_nodes,
                               min_ego_edges = min_ego_edges)

print(results$netdis)
print(results$comp_spec)

## ------------------------------------------------------------------------

binning_fn <- purrr::partial(binned_densities_adaptive,
                             min_counts_per_interval = 10,
                             num_intervals = 50)


# Calculate netdis statistics
results <- netdis_many_to_many(graphs,
                               ref_graph,
                               max_graphlet_size = max_graphlet_size,
                               neighbourhood_size = neighbourhood_size,
                               min_ego_nodes = min_ego_nodes,
                               min_ego_edges = min_ego_edges,
                               binning_fn = binning_fn)

print(results$netdis)
print(results$comp_spec)



## ------------------------------------------------------------------------
bin_counts_fn <- density_binned_counts_gp

exp_counts_fn <- purrr::partial(netdis_expected_counts,
                                scale_fn = NULL)

# Calculate netdis statistics
results <- netdis_many_to_many(graphs,
                               ref_graph = NULL,
                               max_graphlet_size = max_graphlet_size,
                               neighbourhood_size = neighbourhood_size,
                               min_ego_nodes = min_ego_nodes,
                               min_ego_edges = min_ego_edges,
                               bin_counts_fn = bin_counts_fn,
                               exp_counts_fn = exp_counts_fn)

print(results$netdis)
print(results$comp_spec)

## ------------------------------------------------------------------------
binning_fn <- single_density_bin
bin_counts_fn <- density_binned_counts
exp_counts_fn <- netdis_expected_counts

# Calculate netdis statistics
results <- netdis_many_to_many(graphs,
                               ref_graph = NULL,
                               max_graphlet_size = max_graphlet_size,
                               neighbourhood_size = neighbourhood_size,
                               min_ego_nodes = min_ego_nodes,
                               min_ego_edges = min_ego_edges,
                               binning_fn = binning_fn,
                               bin_counts_fn = bin_counts_fn,
                               exp_counts_fn = exp_counts_fn)

print(results$netdis)
print(results$comp_spec)

