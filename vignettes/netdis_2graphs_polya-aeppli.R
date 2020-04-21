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

# Get ego-network densities
densities_1 <- ego_network_density(graphlet_counts_1)
densities_2 <- ego_network_density(graphlet_counts_2)

# Adaptively bin ego-network densities
binned_densities_1 <- binned_densities_adaptive(densities_1, 
                                                min_counts_per_interval = min_bin_count, 
                                                num_intervals = num_bins)

ego_density_bins_1 <- binned_densities_1$breaks

binned_densities_2 <- binned_densities_adaptive(densities_2, 
                                                min_counts_per_interval = min_bin_count, 
                                                num_intervals = num_bins)

ego_density_bins_2 <- binned_densities_2$breaks

## ------------------------------------------------------------------------

#' INTERNAL FUNCTION - DO NOT CALL DIRECTLY
#' Calculate expected counts with geometric poisson (Polya-Aeppli)
#' approximation for a single density bin.
#' @param bin_idx Density bin index to calculate expected counts for.
#' @param graphlet_counts Graphlet counts for a number of ego_networks.
#' @param density_interval_indexes Density bin index for
#' each ego network.
exp_counts_bin_gp <- function(bin_idx, graphlet_counts,
                              density_interval_indexes,
                              mean_binned_graphlet_counts,
                              max_graphlet_size) {
  counts <- graphlet_counts[density_interval_indexes == bin_idx, ]
  means <- mean_binned_graphlet_counts[bin_idx, ]
  
  mean_sub_counts <- sweep(counts, 2, means)
  
  Vd_sq <- colSums(mean_sub_counts^2) / (nrow(mean_sub_counts) - 1)
  theta_d <- 2 * means / (Vd_sq + means)
  
  exp_counts_dk <- vector()
  for (k in 2:max_graphlet_size) {
    graphlet_idx <- graphlet_ids_for_size(k)
    
    lambda_dk <- mean(2 * means[graphlet_idx]^2 /
                        (Vd_sq[graphlet_idx] + means[graphlet_idx]),
                      na.rm = TRUE)
    
    exp_counts_dk <- append(exp_counts_dk,
                            lambda_dk / theta_d[graphlet_idx])
  }
  
  exp_counts_dk
}

#' Calculate expected counts in density bins using the
#' geometric poisson (Polya-Aeppli) approximation.
#' @param graphlet_counts Graphlet counts for a number of ego_networks.
#' @param density_interval_indexes Density bin index for
#' each ego network.
#' @param max_graphlet_size Determines the maximum size of graphlets
#' included in graphlet_counts.
#' @export
density_binned_counts_gp <- function(graphlet_counts,
                                     density_interval_indexes,
                                     max_graphlet_size) {

  mean_binned_graphlet_counts <- mean_density_binned_graphlet_counts(
    graphlet_counts,
    density_interval_indexes)

  nbins <- length(unique(density_interval_indexes))
  expected_counts_bin <- t(sapply(1:nbins,
                                  exp_counts_bin_gp,
                                  graphlet_counts = graphlet_counts,
                                  density_interval_indexes = density_interval_indexes,
                                  mean_binned_graphlet_counts = mean_binned_graphlet_counts,
                                  max_graphlet_size = max_graphlet_size))

  # deal with NAs caused by bins with zero counts for a graphlet
  expected_counts_bin[is.nan(expected_counts_bin)] <- 0

  expected_counts_bin
}

binned_graphlet_counts_1 <- density_binned_counts_gp(graphlet_counts_1,
                                                     binned_densities_1$interval_indexes,
                                                     max_graphlet_size)

binned_graphlet_counts_2 <- density_binned_counts_gp(graphlet_counts_2,
                                                     binned_densities_2$interval_indexes,
                                                     max_graphlet_size)

## ------------------------------------------------------------------------
# Calculate expected graphlet counts for each ego network
exp_graphlet_counts_1 <- netdis_expected_counts(graphlet_counts_1, 
                                                                 ego_density_bins_1, 
                                                                 binned_graphlet_counts_1,
                                                                 max_graphlet_size,
                                                                 scale_fn = NULL)


exp_graphlet_counts_2 <- netdis_expected_counts(graphlet_counts_2, 
                                                                 ego_density_bins_2, 
                                                                 binned_graphlet_counts_2,
                                                                 max_graphlet_size,
                                                                 scale_fn = NULL)
# Centre graphlet counts by subtracting expected counts
centred_graphlet_counts_1 <- netdis_subtract_exp_counts(graphlet_counts_1,
                                                        exp_graphlet_counts_1,
                                                        max_graphlet_size)

centred_graphlet_counts_2 <- netdis_subtract_exp_counts(graphlet_counts_2,
                                                        exp_graphlet_counts_2,
                                                        max_graphlet_size)

## ------------------------------------------------------------------------
sum_graphlet_counts_1 <- colSums(centred_graphlet_counts_1)

sum_graphlet_counts_2 <- colSums(centred_graphlet_counts_2)

## ------------------------------------------------------------------------

netdis_result <- netdis_uptok(sum_graphlet_counts_1, 
                              sum_graphlet_counts_2, 
                              max_graphlet_size)

print(netdis_result)

