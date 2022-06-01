
context("Measures Netdis: Graphlet tuples")
test_message <-
  paste("count_graphlet_tuples and count_graphlet_tuples_ego give",
    "choose(node_count, graphlet_size) for each graph + graphlet",
    "combination",
    sep = " "
  )
test_that(test_message, {
  # Create some test graphs with known node counts (this is the only graph
  # property we care about for this test)
  graph_n11 <- igraph::erdos.renyi.game(11, p = 1, type = "gnp")
  graph_n37 <- igraph::erdos.renyi.game(37, p = 1, type = "gnp")
  graph_n73 <- igraph::erdos.renyi.game(73, p = 1, type = "gnp")

  # calculate graphlet counts object using previously tested function
  graphlet_counts_n11 <- count_graphlets_for_graph(graph_n11, 5)
  graphlet_counts_n37 <- count_graphlets_for_graph(graph_n37, 5)
  graphlet_counts_n73 <- count_graphlets_for_graph(graph_n73, 5)

  # Calculate expected graph tuple count for graphlets of various sizes. There
  # is 1 graphlet of size 1, 2 of size 3, 6 of size 4, and 21 of size 5
  graphlet_tuple_counts <- function(n, max_graphlet_size) {
    if (max_graphlet_size >= 2) {
      tuple_counts <- rep(choose(n, 2), 1)
    }
    if (max_graphlet_size >= 3) {
      tuple_counts <- c(tuple_counts, rep(choose(n, 3), 2))
    }
    if (max_graphlet_size >= 4) {
      tuple_counts <- c(tuple_counts, rep(choose(n, 4), 6))
    }
    if (max_graphlet_size >= 5) {
      tuple_counts <- c(tuple_counts, rep(choose(n, 5), 21))
    }
    tuple_counts <- setNames(tuple_counts, graphlet_key(max_graphlet_size)$id)
    tuple_counts
  }

  # === TEST count_graphlet_tuples ===
  # Generate expected tuple counts for graphlets up to size 4 and 5
  expected_tuple_count_n11_gs4 <- graphlet_tuple_counts(11, 4)
  expected_tuple_count_n37_gs4 <- graphlet_tuple_counts(37, 4)
  expected_tuple_count_n73_gs4 <- graphlet_tuple_counts(73, 4)
  expected_tuple_count_n11_gs5 <- graphlet_tuple_counts(11, 5)
  expected_tuple_count_n37_gs5 <- graphlet_tuple_counts(37, 5)
  expected_tuple_count_n73_gs5 <- graphlet_tuple_counts(73, 5)

  # Generate actual tuple counts for graphlets up to size 4 and 5
  actual_tuple_count_n11_gs4 <- count_graphlet_tuples(graphlet_counts_n11, 4)
  actual_tuple_count_n37_gs4 <- count_graphlet_tuples(graphlet_counts_n37, 4)
  actual_tuple_count_n73_gs4 <- count_graphlet_tuples(graphlet_counts_n73, 4)
  actual_tuple_count_n11_gs5 <- count_graphlet_tuples(graphlet_counts_n11, 5)
  actual_tuple_count_n37_gs5 <- count_graphlet_tuples(graphlet_counts_n37, 5)
  actual_tuple_count_n73_gs5 <- count_graphlet_tuples(graphlet_counts_n73, 5)

  # Compare expected tuple counts with actual
  expect_equal(expected_tuple_count_n11_gs4, actual_tuple_count_n11_gs4)
  expect_equal(expected_tuple_count_n37_gs4, actual_tuple_count_n37_gs4)
  expect_equal(expected_tuple_count_n73_gs4, actual_tuple_count_n73_gs4)
  expect_equal(expected_tuple_count_n11_gs5, actual_tuple_count_n11_gs5)
  expect_equal(expected_tuple_count_n37_gs5, actual_tuple_count_n37_gs5)
  expect_equal(expected_tuple_count_n73_gs5, actual_tuple_count_n73_gs5)

  # === TEST count_graphlet_tuples_ego ===
  # NOTE: This test is not amazing, as graphlet_tuple_counts_ego is very similar
  # to the method under test. However, it's a simple method so maybe that's ok?
  graphlet_tuple_counts_ego <- function(graphlet_counts_ego, max_graphlet_size) {
    t(apply(graphlet_counts_ego, 1,
      count_graphlet_tuples,
      max_graphlet_size = max_graphlet_size
    ))
  }
  # Generate ego network graphlet counts for each graph
  graph_n11_ego1 <- count_graphlets_ego(graph_n11, neighbourhood_size = 1)
  graph_n37_ego1 <- count_graphlets_ego(graph_n37, neighbourhood_size = 1)
  graph_n73_ego1 <- count_graphlets_ego(graph_n73, neighbourhood_size = 1)
  graph_n11_ego2 <- count_graphlets_ego(graph_n11, neighbourhood_size = 2)
  graph_n37_ego2 <- count_graphlets_ego(graph_n37, neighbourhood_size = 2)
  graph_n73_ego2 <- count_graphlets_ego(graph_n73, neighbourhood_size = 2)

  # Generate expected tuple counts for graphlets up to size 4 and 5
  # 1. For ego-networks of order 1
  expected_tuple_count_n11_ego1_gs4 <- graphlet_tuple_counts_ego(graph_n11_ego1, 4)
  expected_tuple_count_n37_ego1_gs4 <- graphlet_tuple_counts_ego(graph_n37_ego1, 4)
  expected_tuple_count_n73_ego1_gs4 <- graphlet_tuple_counts_ego(graph_n73_ego1, 4)
  expected_tuple_count_n11_ego1_gs5 <- graphlet_tuple_counts_ego(graph_n11_ego1, 5)
  expected_tuple_count_n37_ego1_gs5 <- graphlet_tuple_counts_ego(graph_n37_ego1, 5)
  expected_tuple_count_n73_ego1_gs5 <- graphlet_tuple_counts_ego(graph_n73_ego1, 5)
  # 2. For ego-networks of order 2
  expected_tuple_count_n11_ego2_gs4 <- graphlet_tuple_counts_ego(graph_n11_ego2, 4)
  expected_tuple_count_n37_ego2_gs4 <- graphlet_tuple_counts_ego(graph_n37_ego2, 4)
  expected_tuple_count_n73_ego2_gs4 <- graphlet_tuple_counts_ego(graph_n73_ego2, 4)
  expected_tuple_count_n11_ego2_gs5 <- graphlet_tuple_counts_ego(graph_n11_ego2, 5)
  expected_tuple_count_n37_ego2_gs5 <- graphlet_tuple_counts_ego(graph_n37_ego2, 5)
  expected_tuple_count_n73_ego2_gs5 <- graphlet_tuple_counts_ego(graph_n73_ego2, 5)

  # Calculate actual tuple counts
  # 1. For ego-networks of order 1
  actual_tuple_count_n11_ego1_gs4 <- count_graphlet_tuples_ego(graph_n11_ego1, 4)
  actual_tuple_count_n37_ego1_gs4 <- count_graphlet_tuples_ego(graph_n37_ego1, 4)
  actual_tuple_count_n73_ego1_gs4 <- count_graphlet_tuples_ego(graph_n73_ego1, 4)
  actual_tuple_count_n11_ego1_gs5 <- count_graphlet_tuples_ego(graph_n11_ego1, 5)
  actual_tuple_count_n37_ego1_gs5 <- count_graphlet_tuples_ego(graph_n37_ego1, 5)
  actual_tuple_count_n73_ego1_gs5 <- count_graphlet_tuples_ego(graph_n73_ego1, 5)
  # 2. For ego-networks of order 2
  actual_tuple_count_n11_ego2_gs4 <- count_graphlet_tuples_ego(graph_n11_ego2, 4)
  actual_tuple_count_n37_ego2_gs4 <- count_graphlet_tuples_ego(graph_n37_ego2, 4)
  actual_tuple_count_n73_ego2_gs4 <- count_graphlet_tuples_ego(graph_n73_ego2, 4)
  actual_tuple_count_n11_ego2_gs5 <- count_graphlet_tuples_ego(graph_n11_ego2, 5)
  actual_tuple_count_n37_ego2_gs5 <- count_graphlet_tuples_ego(graph_n37_ego2, 5)
  actual_tuple_count_n73_ego2_gs5 <- count_graphlet_tuples_ego(graph_n73_ego2, 5)

  # Compare expected with actual
  expect_equal(expected_tuple_count_n11_ego1_gs4, actual_tuple_count_n11_ego1_gs4)
  expect_equal(expected_tuple_count_n37_ego1_gs4, actual_tuple_count_n37_ego1_gs4)
  expect_equal(expected_tuple_count_n73_ego1_gs4, actual_tuple_count_n73_ego1_gs4)
  expect_equal(expected_tuple_count_n11_ego1_gs5, actual_tuple_count_n11_ego1_gs5)
  expect_equal(expected_tuple_count_n37_ego1_gs5, actual_tuple_count_n37_ego1_gs5)
  expect_equal(expected_tuple_count_n73_ego1_gs5, actual_tuple_count_n73_ego1_gs5)
  # 2. For ego-networks of order 2
  expect_equal(expected_tuple_count_n11_ego2_gs4, actual_tuple_count_n11_ego2_gs4)
  expect_equal(expected_tuple_count_n37_ego2_gs4, actual_tuple_count_n37_ego2_gs4)
  expect_equal(expected_tuple_count_n73_ego2_gs4, actual_tuple_count_n73_ego2_gs4)
  expect_equal(expected_tuple_count_n11_ego2_gs5, actual_tuple_count_n11_ego2_gs5)
  expect_equal(expected_tuple_count_n37_ego2_gs5, actual_tuple_count_n37_ego2_gs5)
  expect_equal(expected_tuple_count_n73_ego2_gs5, actual_tuple_count_n73_ego2_gs5)
})

context("Measures Netdis: Ego-network density values match those for manually verified networks")
test_that("Ego-network 4-node density values match manually verified totals", {
  # Set up a small sample network with at least one ego-network that contains
  # at least one of each graphlets
  elist <- rbind(
    c("n1", "n2"),
    c("n2", "n3"),
    c("n1", "n4"),
    c("n2", "n5"),
    c("n1", "n6"),
    c("n1", "n7"),
    c("n2", "n4"),
    c("n4", "n6"),
    c("n6", "n8"),
    c("n7", "n8"),
    c("n7", "n9"),
    c("n7", "n10"),
    c("n8", "n9"),
    c("n8", "n10"),
    c("n9", "n10")
  )
  graph <- igraph::graph_from_edgelist(elist, directed = FALSE)

  # Set parameters for test
  max_graphlet_size <- 4
  min_counts_per_interval <- 2
  num_intervals <- 100
  min_ego_edges <- 0
  min_ego_nodes <- 0

  # Set node and graphlet labels to use for row and col names in expected counts
  node_labels <- igraph::V(graph)$name
  graphlet_labels <- c("G0", "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8")

  # Set manually verified ego-network node counts and edge densities
  # 1 . Ego-networks of order 1
  expected_node_counts_o1 <- c(5, 5, 2, 4, 2, 4, 5, 5, 4, 4)
  expected_edge_counts_o1 <- c(6, 5, 1, 5, 1, 4, 7, 7, 6, 6)
  max_edge_counts_o1 <- choose(expected_node_counts_o1, 2)
  expected_densities_o1 <- c(expected_edge_counts_o1 / max_edge_counts_o1)
  names(expected_densities_o1) <- node_labels
  # Order 1 expected densities should be:
  # 0.6, 0.5, 1.0, 0.83, 1.0, 0.67, 0.7, 0.7, 1.0, 1.0
  # 2. Ego-networks of order 2
  expected_node_counts_o2 <- c(10, 7, 5, 8, 5, 8, 8, 7, 6, 6)
  expected_edge_counts_o2 <- c(15, 8, 5, 10, 5, 13, 13, 11, 9, 9)
  max_edge_counts_o2 <- choose(expected_node_counts_o2, 2)
  expected_densities_o2 <- c(expected_edge_counts_o2 / max_edge_counts_o2)
  names(expected_densities_o2) <- node_labels
  # Order 2 expected densities should be:
  # 0.33, 0.38, 0.50, 0.36, 0.50, 0.46, 0.46, 0.52, 0.60, 0.60

  # Generate order 1 and 2 ego network graphlet counts
  # with previously tested function
  graphlet_counts_ego_o1 <- count_graphlets_ego(graph,
    neighbourhood_size = 1,
    min_ego_edges = min_ego_edges,
    min_ego_nodes = min_ego_nodes
  )
  graphlet_counts_ego_o2 <- count_graphlets_ego(graph,
    neighbourhood_size = 2,
    min_ego_edges = min_ego_edges,
    min_ego_nodes = min_ego_nodes
  )

  # Calculate densities
  actual_densities_o1 <- ego_network_density(graphlet_counts_ego_o1)
  actual_densities_o2 <- ego_network_density(graphlet_counts_ego_o2)

  # Check densities match expected values
  expect_equal(actual_densities_o1, expected_densities_o1)
  expect_equal(actual_densities_o2, expected_densities_o2)
})

context("Measures Netdis: Ego-network density-binned reference counts for manually verified networks")
test_that("Ego-network 4-node density-binned reference counts match manually verified totals", {
  # Set up a small sample network with at least one ego-network that contains
  # at least one of each graphlets
  elist <- rbind(
    c("n1", "n2"),
    c("n2", "n3"),
    c("n1", "n4"),
    c("n2", "n5"),
    c("n1", "n6"),
    c("n1", "n7"),
    c("n2", "n4"),
    c("n4", "n6"),
    c("n6", "n8"),
    c("n7", "n8"),
    c("n7", "n9"),
    c("n7", "n10"),
    c("n8", "n9"),
    c("n8", "n10"),
    c("n9", "n10")
  )
  graph <- igraph::graph_from_edgelist(elist, directed = FALSE)

  # Set parameters for test
  max_graphlet_size <- 4
  min_counts_per_interval <- 2
  num_intervals <- 100

  # Set node and graphlet labels to use for row and col names in expected counts
  node_labels <- igraph::V(graph)$name
  graphlet_labels <- c("G0", "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8")

  # Set manually verified ego-network node counts and edge densities
  # 1 . Ego-networks of order 1
  expected_node_counts_o1 <- c(5, 5, 2, 4, 2, 4, 5, 5, 4, 4)
  expected_edge_counts_o1 <- c(6, 5, 1, 5, 1, 4, 7, 7, 6, 6)
  max_edge_counts_o1 <- choose(expected_node_counts_o1, 2)
  expected_densities_o1 <- c(expected_edge_counts_o1 / max_edge_counts_o1)
  # Order 1 expected densities should be:
  # 0.6, 0.5, 1.0, 0.83, 1.0, 0.67, 0.7, 0.7, 1.0, 1.0
  # 2. Ego-networks of order 2
  expected_node_counts_o2 <- c(10, 7, 5, 8, 5, 8, 8, 7, 6, 6)
  expected_edge_counts_o2 <- c(15, 8, 5, 10, 5, 13, 13, 11, 9, 9)
  max_edge_counts_o2 <- choose(expected_node_counts_o2, 2)
  expected_densities_o2 <- c(expected_edge_counts_o2 / max_edge_counts_o2)
  # Order 2 expected densities should be:
  # 0.33, 0.38, 0.50, 0.36, 0.50, 0.46, 0.46, 0.52, 0.60, 0.60

  # Set manually verified density bins for ego-networks
  # 1. Ego-networks of order 1
  expected_breaks_o1 <- c(0.5, 0.605, 0.705, 1)
  expected_interval_indexes_o1 <- c(1, 1, 3, 3, 3, 2, 2, 2, 3, 3)
  expected_binned_densities_o1 <- list(
    densities = expected_densities_o1,
    interval_indexes = expected_interval_indexes_o1,
    breaks = expected_breaks_o1
  )
  # Check binned densities are as expected
  actual_binned_densities_o1 <- binned_densities_adaptive(
    expected_densities_o1,
    min_counts_per_interval = min_counts_per_interval,
    num_intervals = num_intervals
  )
  expect_equal(actual_binned_densities_o1, expected_binned_densities_o1)
  # 2. Ego-networks of order 2
  expected_min_break_o2 <- 1 / 3
  expected_max_break_o2 <- 0.6
  expected_initial_interval_o2 <-
    (expected_max_break_o2 - expected_min_break_o2) / (num_intervals) # 0.00266666667
  expected_breaks_o2 <- expected_min_break_o2 + (expected_initial_interval_o2 * c(0, 9, 50, 63, 100))
  expected_interval_indexes_o2 <- c(1, 2, 3, 1, 3, 2, 2, 4, 4, 4)
  expected_binned_densities_o2 <- list(
    densities = expected_densities_o2,
    interval_indexes = expected_interval_indexes_o2,
    breaks = expected_breaks_o2
  )
  # Check binned densities are as expected
  actual_binned_densities_o2 <- binned_densities_adaptive(
    expected_densities_o2,
    min_counts_per_interval = min_counts_per_interval,
    num_intervals = num_intervals
  )
  expect_equal(actual_binned_densities_o2, expected_binned_densities_o2)

  # Set manually verified scaled ego-network graphlet counts
  graphlet_key <- graphlet_key(max_graphlet_size)
  k <- graphlet_key$node_count
  # 1-step ego networks
  expected_counts_o1 <- rbind(
    c(6, 5, 2, 0, 1, 0, 2, 1, 0) / zeros_to_ones(choose(5, k)),
    c(5, 5, 1, 0, 2, 0, 2, 0, 0) / zeros_to_ones(choose(5, k)),
    c(1, 0, 0, 0, 0, 0, 0, 0, 0) / zeros_to_ones(choose(2, k)),
    c(5, 2, 2, 0, 0, 0, 0, 1, 0) / zeros_to_ones(choose(4, k)),
    c(1, 0, 0, 0, 0, 0, 0, 0, 0) / zeros_to_ones(choose(2, k)),
    c(4, 2, 1, 0, 0, 0, 1, 0, 0) / zeros_to_ones(choose(4, k)),
    c(7, 3, 4, 0, 0, 0, 3, 0, 1) / zeros_to_ones(choose(5, k)),
    c(7, 3, 4, 0, 0, 0, 3, 0, 1) / zeros_to_ones(choose(5, k)),
    c(6, 0, 4, 0, 0, 0, 0, 0, 1) / zeros_to_ones(choose(4, k)),
    c(6, 0, 4, 0, 0, 0, 0, 0, 1) / zeros_to_ones(choose(4, k))
  )
  rownames(expected_counts_o1) <- node_labels
  colnames(expected_counts_o1) <- graphlet_labels
  # 2-step ego networks
  expected_counts_o2 <- rbind(
    c(15, 18, 6, 21, 3, 1, 11, 1, 1) / zeros_to_ones(choose(10, k)),
    c(8, 10, 2, 6, 3, 0, 4, 1, 0) / zeros_to_ones(choose(7, k)),
    c(5, 5, 1, 0, 2, 0, 2, 0, 0) / zeros_to_ones(choose(5, k)),
    c(10, 14, 2, 11, 3, 1, 5, 1, 0) / zeros_to_ones(choose(8, k)),
    c(5, 5, 1, 0, 2, 0, 2, 0, 0) / zeros_to_ones(choose(5, k)),
    c(13, 13, 6, 15, 1, 1, 9, 1, 1) / zeros_to_ones(choose(8, k)),
    c(13, 13, 6, 15, 1, 1, 9, 1, 1) / zeros_to_ones(choose(8, k)),
    c(11, 10, 5, 10, 0, 1, 8, 0, 1) / zeros_to_ones(choose(7, k)),
    c(9, 8, 4, 4, 0, 1, 6, 0, 1) / zeros_to_ones(choose(6, k)),
    c(9, 8, 4, 4, 0, 1, 6, 0, 1) / zeros_to_ones(choose(6, k))
  )
  rownames(expected_counts_o2) <- node_labels
  colnames(expected_counts_o2) <- graphlet_labels

  # Calculate binned average expected counts based on manually verified counts
  # and density bins
  # Order 1: Expected interval indexes: 1, 1, 3, 3, 3, 2, 2, 2, 3, 3
  mean_counts_bin1_o1 <- (expected_counts_o1[1, ] + expected_counts_o1[2, ]) / 2
  mean_counts_bin2_o1 <- (expected_counts_o1[6, ] + expected_counts_o1[7, ] +
    expected_counts_o1[8, ]) / 3
  mean_counts_bin3_o1 <- (expected_counts_o1[3, ] + expected_counts_o1[4, ] +
    expected_counts_o1[5, ] + expected_counts_o1[9, ] +
    expected_counts_o1[10, ]) / 5
  expected_mean_density_binned_counts_o1 <- rbind(
    mean_counts_bin1_o1, mean_counts_bin2_o1, mean_counts_bin3_o1
  )
  rownames(expected_mean_density_binned_counts_o1) <- 1:3
  # Order 2: Expected interval indexes: 1, 3, 3, 1, 3, 2, 2, 4, 4, 4
  mean_counts_bin1_o2 <- (expected_counts_o2[1, ] + expected_counts_o2[4, ]) / 2
  mean_counts_bin2_o2 <- (expected_counts_o2[2, ] + expected_counts_o2[6, ] +
    expected_counts_o2[7, ]) / 3
  mean_counts_bin3_o2 <- (expected_counts_o2[3, ] + expected_counts_o2[5, ]) / 2
  mean_counts_bin4_o2 <- (expected_counts_o2[8, ] + expected_counts_o2[9, ] +
    expected_counts_o2[10, ]) / 3
  expected_mean_density_binned_counts_o2 <- rbind(
    mean_counts_bin1_o2, mean_counts_bin2_o2, mean_counts_bin3_o2,
    mean_counts_bin4_o2
  )
  rownames(expected_mean_density_binned_counts_o2) <- 1:4

  # Calculate actual output of function under test
  actual_mean_density_binned_counts_o1 <- mean_density_binned_graphlet_counts(
    expected_counts_o1, expected_interval_indexes_o1
  )
  actual_mean_density_binned_counts_o2 <- mean_density_binned_graphlet_counts(
    expected_counts_o2, expected_interval_indexes_o2
  )

  # Check actual output vs expected
  expect_equal(
    actual_mean_density_binned_counts_o1,
    expected_mean_density_binned_counts_o1
  )
  expect_equal(
    actual_mean_density_binned_counts_o2,
    expected_mean_density_binned_counts_o2
  )
})

context("Measures Netdis: scale_graphlet_counts_ego for manually verified networks")
test_that("Ego-network 4-node graphlet counts match manually verified totals", {
  # Set up a small sample network with at least one ego-network that contains
  # at least one of each graphlets
  elist <- rbind(
    c("n1", "n2"),
    c("n2", "n3"),
    c("n1", "n4"),
    c("n2", "n5"),
    c("n1", "n6"),
    c("n1", "n7"),
    c("n2", "n4"),
    c("n4", "n6"),
    c("n6", "n8"),
    c("n7", "n8"),
    c("n7", "n9"),
    c("n7", "n10"),
    c("n8", "n9"),
    c("n8", "n10"),
    c("n9", "n10")
  )
  graph <- igraph::graph_from_edgelist(elist, directed = FALSE)

  # Set node and graphlet labels to use for row and col names in expected counts
  node_labels <- igraph::V(graph)$name
  graphlet_labels <- c("G0", "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8")

  # Count graphlets in each ego network of the graph with neighbourhood sizes of 1 and 2
  max_graphlet_size <- 4
  min_ego_edges <- 0
  min_ego_nodes <- 0

  # Use previously tested functions to generate ego networks and calcualte graphlet
  # counts.
  # ego nets
  ego_networks_o1 <- make_named_ego_graph(graph,
    order = 1,
    min_ego_edges = min_ego_edges,
    min_ego_nodes = min_ego_nodes
  )
  ego_networks_o2 <- make_named_ego_graph(graph,
    order = 2,
    min_ego_edges = min_ego_edges,
    min_ego_nodes = min_ego_nodes
  )

  # graphlet counts
  graphlet_counts_o1 <-
    ego_to_graphlet_counts(ego_networks_o1,
      max_graphlet_size = max_graphlet_size
    )
  graphlet_counts_o2 <-
    ego_to_graphlet_counts(ego_networks_o2,
      max_graphlet_size = max_graphlet_size
    )


  # Calculate scaled counts with scale_graphlet_counts_ego
  # (function to test).
  actual_counts_o1 <-
    scale_graphlet_counts_ego(graphlet_counts_o1,
      max_graphlet_size = max_graphlet_size
    )
  actual_counts_o2 <-
    scale_graphlet_counts_ego(graphlet_counts_o2,
      max_graphlet_size = max_graphlet_size
    )

  graphlet_key <- graphlet_key(max_graphlet_size)
  k <- graphlet_key$node_count
  # Set manually verified counts
  # 1-step ego networks
  expected_counts_o1 <- rbind(
    c(6, 5, 2, 0, 1, 0, 2, 1, 0) / zeros_to_ones(choose(5, k)),
    c(5, 5, 1, 0, 2, 0, 2, 0, 0) / zeros_to_ones(choose(5, k)),
    c(1, 0, 0, 0, 0, 0, 0, 0, 0) / zeros_to_ones(choose(2, k)),
    c(5, 2, 2, 0, 0, 0, 0, 1, 0) / zeros_to_ones(choose(4, k)),
    c(1, 0, 0, 0, 0, 0, 0, 0, 0) / zeros_to_ones(choose(2, k)),
    c(4, 2, 1, 0, 0, 0, 1, 0, 0) / zeros_to_ones(choose(4, k)),
    c(7, 3, 4, 0, 0, 0, 3, 0, 1) / zeros_to_ones(choose(5, k)),
    c(7, 3, 4, 0, 0, 0, 3, 0, 1) / zeros_to_ones(choose(5, k)),
    c(6, 0, 4, 0, 0, 0, 0, 0, 1) / zeros_to_ones(choose(4, k)),
    c(6, 0, 4, 0, 0, 0, 0, 0, 1) / zeros_to_ones(choose(4, k))
  )
  rownames(expected_counts_o1) <- node_labels
  colnames(expected_counts_o1) <- graphlet_labels
  # 2-step ego networks
  expected_counts_o2 <- rbind(
    c(15, 18, 6, 21, 3, 1, 11, 1, 1) / zeros_to_ones(choose(10, k)),
    c(8, 10, 2, 6, 3, 0, 4, 1, 0) / zeros_to_ones(choose(7, k)),
    c(5, 5, 1, 0, 2, 0, 2, 0, 0) / zeros_to_ones(choose(5, k)),
    c(10, 14, 2, 11, 3, 1, 5, 1, 0) / zeros_to_ones(choose(8, k)),
    c(5, 5, 1, 0, 2, 0, 2, 0, 0) / zeros_to_ones(choose(5, k)),
    c(13, 13, 6, 15, 1, 1, 9, 1, 1) / zeros_to_ones(choose(8, k)),
    c(13, 13, 6, 15, 1, 1, 9, 1, 1) / zeros_to_ones(choose(8, k)),
    c(11, 10, 5, 10, 0, 1, 8, 0, 1) / zeros_to_ones(choose(7, k)),
    c(9, 8, 4, 4, 0, 1, 6, 0, 1) / zeros_to_ones(choose(6, k)),
    c(9, 8, 4, 4, 0, 1, 6, 0, 1) / zeros_to_ones(choose(6, k))
  )
  rownames(expected_counts_o2) <- node_labels
  colnames(expected_counts_o2) <- graphlet_labels

  # Test that actual counts match expected
  expect_equal(actual_counts_o1, expected_counts_o1)
  expect_equal(actual_counts_o2, expected_counts_o2)
})

context("Measures Netdis: Ego-network density-binned counts for manually verified networks")
test_that("density_binned_counts output matches manually verified totals with different scaling and aggregation functions", {
  # Set up a small sample network with at least one ego-network that contains
  # at least one of each graphlets
  elist <- rbind(
    c("n1", "n2"),
    c("n2", "n3"),
    c("n1", "n4"),
    c("n2", "n5"),
    c("n1", "n6"),
    c("n1", "n7"),
    c("n2", "n4"),
    c("n4", "n6"),
    c("n6", "n8"),
    c("n7", "n8"),
    c("n7", "n9"),
    c("n7", "n10"),
    c("n8", "n9"),
    c("n8", "n10"),
    c("n9", "n10")
  )
  graph <- igraph::graph_from_edgelist(elist, directed = FALSE)

  # Set parameters for test
  max_graphlet_size <- 4
  min_counts_per_interval <- 2
  num_intervals <- 100
  min_ego_edges <- 0
  min_ego_nodes <- 0

  # Set node and graphlet labels to use for row and col names in expected counts
  node_labels <- igraph::V(graph)$name
  graphlet_labels <- c("G0", "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8")

  # Set manually verified ego-network node counts and edge densities
  # 1 . Ego-networks of order 1
  expected_node_counts_o1 <- c(5, 5, 2, 4, 2, 4, 5, 5, 4, 4)
  expected_edge_counts_o1 <- c(6, 5, 1, 5, 1, 4, 7, 7, 6, 6)
  max_edge_counts_o1 <- choose(expected_node_counts_o1, 2)
  expected_densities_o1 <- c(expected_edge_counts_o1 / max_edge_counts_o1)
  # Order 1 expected densities should be:
  # 0.6, 0.5, 1.0, 0.83, 1.0, 0.67, 0.7, 0.7, 1.0, 1.0
  # 2. Ego-networks of order 2
  expected_node_counts_o2 <- c(10, 7, 5, 8, 5, 8, 8, 7, 6, 6)
  expected_edge_counts_o2 <- c(15, 8, 5, 10, 5, 13, 13, 11, 9, 9)
  max_edge_counts_o2 <- choose(expected_node_counts_o2, 2)
  expected_densities_o2 <- c(expected_edge_counts_o2 / max_edge_counts_o2)
  # Order 2 expected densities should be:
  # 0.33, 0.38, 0.50, 0.36, 0.50, 0.46, 0.46, 0.52, 0.60, 0.60

  # Set manually verified density bins for ego-networks
  # 1. Ego-networks of order 1
  expected_breaks_o1 <- c(0.5, 0.605, 0.705, 1)
  expected_interval_indexes_o1 <- c(1, 1, 3, 3, 3, 2, 2, 2, 3, 3)
  expected_binned_densities_o1 <- list(
    densities = expected_densities_o1,
    interval_indexes = expected_interval_indexes_o1,
    breaks = expected_breaks_o1
  )
  # 2. Ego-networks of order 2
  expected_min_break_o2 <- 1 / 3
  expected_max_break_o2 <- 0.6
  expected_initial_interval_o2 <-
    (expected_max_break_o2 - expected_min_break_o2) / (num_intervals) # 0.00266666667
  expected_breaks_o2 <- expected_min_break_o2 + (expected_initial_interval_o2 * c(0, 9, 50, 63, 100))
  expected_interval_indexes_o2 <- c(1, 2, 3, 1, 3, 2, 2, 4, 4, 4)
  expected_binned_densities_o2 <- list(
    densities = expected_densities_o2,
    interval_indexes = expected_interval_indexes_o2,
    breaks = expected_breaks_o2
  )

  # Set manually verified scaled ego-network graphlet counts
  graphlet_key <- graphlet_key(max_graphlet_size)
  k <- graphlet_key$node_count
  # 1-step ego networks
  expected_counts_o1 <- rbind(
    c(6, 5, 2, 0, 1, 0, 2, 1, 0) / zeros_to_ones(choose(5, k)),
    c(5, 5, 1, 0, 2, 0, 2, 0, 0) / zeros_to_ones(choose(5, k)),
    c(1, 0, 0, 0, 0, 0, 0, 0, 0) / zeros_to_ones(choose(2, k)),
    c(5, 2, 2, 0, 0, 0, 0, 1, 0) / zeros_to_ones(choose(4, k)),
    c(1, 0, 0, 0, 0, 0, 0, 0, 0) / zeros_to_ones(choose(2, k)),
    c(4, 2, 1, 0, 0, 0, 1, 0, 0) / zeros_to_ones(choose(4, k)),
    c(7, 3, 4, 0, 0, 0, 3, 0, 1) / zeros_to_ones(choose(5, k)),
    c(7, 3, 4, 0, 0, 0, 3, 0, 1) / zeros_to_ones(choose(5, k)),
    c(6, 0, 4, 0, 0, 0, 0, 0, 1) / zeros_to_ones(choose(4, k)),
    c(6, 0, 4, 0, 0, 0, 0, 0, 1) / zeros_to_ones(choose(4, k))
  )
  rownames(expected_counts_o1) <- node_labels
  colnames(expected_counts_o1) <- graphlet_labels
  # 2-step ego networks
  expected_counts_o2 <- rbind(
    c(15, 18, 6, 21, 3, 1, 11, 1, 1) / zeros_to_ones(choose(10, k)),
    c(8, 10, 2, 6, 3, 0, 4, 1, 0) / zeros_to_ones(choose(7, k)),
    c(5, 5, 1, 0, 2, 0, 2, 0, 0) / zeros_to_ones(choose(5, k)),
    c(10, 14, 2, 11, 3, 1, 5, 1, 0) / zeros_to_ones(choose(8, k)),
    c(5, 5, 1, 0, 2, 0, 2, 0, 0) / zeros_to_ones(choose(5, k)),
    c(13, 13, 6, 15, 1, 1, 9, 1, 1) / zeros_to_ones(choose(8, k)),
    c(13, 13, 6, 15, 1, 1, 9, 1, 1) / zeros_to_ones(choose(8, k)),
    c(11, 10, 5, 10, 0, 1, 8, 0, 1) / zeros_to_ones(choose(7, k)),
    c(9, 8, 4, 4, 0, 1, 6, 0, 1) / zeros_to_ones(choose(6, k)),
    c(9, 8, 4, 4, 0, 1, 6, 0, 1) / zeros_to_ones(choose(6, k))
  )
  rownames(expected_counts_o2) <- node_labels
  colnames(expected_counts_o2) <- graphlet_labels

  # Calculate binned average expected counts based on manually verified counts
  # and density bins
  # Order 1: Expected interval indexes: 1, 1, 3, 3, 3, 2, 2, 2, 3, 3
  mean_counts_bin1_o1 <- (expected_counts_o1[1, ] + expected_counts_o1[2, ]) / 2
  mean_counts_bin2_o1 <- (expected_counts_o1[6, ] + expected_counts_o1[7, ] +
    expected_counts_o1[8, ]) / 3
  mean_counts_bin3_o1 <- (expected_counts_o1[3, ] + expected_counts_o1[4, ] +
    expected_counts_o1[5, ] + expected_counts_o1[9, ] +
    expected_counts_o1[10, ]) / 5
  expected_mean_density_binned_counts_o1 <- rbind(
    mean_counts_bin1_o1, mean_counts_bin2_o1, mean_counts_bin3_o1
  )
  rownames(expected_mean_density_binned_counts_o1) <- 1:3
  # Order 2: Expected interval indexes: 1, 3, 3, 1, 3, 2, 2, 4, 4, 4
  mean_counts_bin1_o2 <- (expected_counts_o2[1, ] + expected_counts_o2[4, ]) / 2
  mean_counts_bin2_o2 <- (expected_counts_o2[2, ] + expected_counts_o2[6, ] +
    expected_counts_o2[7, ]) / 3
  mean_counts_bin3_o2 <- (expected_counts_o2[3, ] + expected_counts_o2[5, ]) / 2
  mean_counts_bin4_o2 <- (expected_counts_o2[8, ] + expected_counts_o2[9, ] +
    expected_counts_o2[10, ]) / 3
  expected_mean_density_binned_counts_o2 <- rbind(
    mean_counts_bin1_o2, mean_counts_bin2_o2, mean_counts_bin3_o2,
    mean_counts_bin4_o2
  )
  rownames(expected_mean_density_binned_counts_o2) <- 1:4

  # density_binned_counts with default arguments should give
  # mean graphlet count in each density bin
  actual_density_binned_counts_o1 <- density_binned_counts(
    expected_counts_o1,
    expected_interval_indexes_o1
  )

  actual_density_binned_counts_o2 <- density_binned_counts(
    expected_counts_o2,
    expected_interval_indexes_o2
  )

  # Check actual output vs expected
  expect_equal(
    actual_density_binned_counts_o1,
    expected_mean_density_binned_counts_o1
  )
  expect_equal(
    actual_density_binned_counts_o2,
    expected_mean_density_binned_counts_o2
  )

  # Calculate max binned counts based on manually verified counts
  # and density bins
  # Order 1: Expected interval indexes: 1, 1, 3, 3, 3, 2, 2, 2, 3, 3
  # apply(x, 2, max): returns max of each column in x
  max_counts_bin1_o1 <- apply(rbind(expected_counts_o1[1, ], expected_counts_o1[2, ]), 2, max)
  max_counts_bin2_o1 <- apply(rbind(
    expected_counts_o1[6, ], expected_counts_o1[7, ],
    expected_counts_o1[8, ]
  ), 2, max)
  max_counts_bin3_o1 <- apply(rbind(
    expected_counts_o1[3, ], expected_counts_o1[4, ],
    expected_counts_o1[5, ], expected_counts_o1[9, ],
    expected_counts_o1[10, ]
  ), 2, max)

  expected_max_density_binned_counts_o1 <- rbind(
    max_counts_bin1_o1, max_counts_bin2_o1, max_counts_bin3_o1
  )
  rownames(expected_max_density_binned_counts_o1) <- 1:3
  # Order 2: Expected interval indexes: 1, 3, 3, 1, 3, 2, 2, 4, 4, 4
  max_counts_bin1_o2 <- apply(rbind(expected_counts_o2[1, ], expected_counts_o2[4, ]), 2, max)
  max_counts_bin2_o2 <- apply(rbind(
    expected_counts_o2[2, ], expected_counts_o2[6, ],
    expected_counts_o2[7, ]
  ), 2, max)
  max_counts_bin3_o2 <- apply(rbind(expected_counts_o2[3, ], expected_counts_o2[5, ]), 2, max)
  max_counts_bin4_o2 <- apply(rbind(
    expected_counts_o2[8, ], expected_counts_o2[9, ],
    expected_counts_o2[10, ]
  ), 2, max)

  expected_max_density_binned_counts_o2 <- rbind(
    max_counts_bin1_o2, max_counts_bin2_o2, max_counts_bin3_o2,
    max_counts_bin4_o2
  )
  rownames(expected_max_density_binned_counts_o2) <- 1:4

  # density_binned_counts with agg_fn = max should give
  # max graphlet count in each density bin
  agg_fn <- max
  scale_fn <- NULL

  actual_max_density_binned_counts_o1 <- density_binned_counts(
    expected_counts_o1,
    expected_interval_indexes_o1,
    agg_fn = agg_fn,
    scale_fn = scale_fn
  )

  actual_max_density_binned_counts_o2 <- density_binned_counts(
    expected_counts_o2,
    expected_interval_indexes_o2,
    agg_fn = agg_fn,
    scale_fn = scale_fn
  )

  # Check actual output vs expected
  expect_equal(
    actual_max_density_binned_counts_o1,
    expected_max_density_binned_counts_o1
  )
  expect_equal(
    actual_max_density_binned_counts_o2,
    expected_max_density_binned_counts_o2
  )

  # density_binned_counts with scale_fn = scale_graphlet_counts_ego
  # should give mean graphlet counts in each density bin scaled by
  # count_graphlet_tuples.
  agg_fn <- mean
  scale_fn <- scale_graphlet_counts_ego

  # calculate expected counts using previously tested function
  expected_scaled_counts_o1 <-
    scale_graphlet_counts_ego(expected_counts_o1,
      max_graphlet_size = max_graphlet_size
    )
  expected_scaled_counts_o2 <-
    scale_graphlet_counts_ego(expected_counts_o2,
      max_graphlet_size = max_graphlet_size
    )

  # calculate mean expected counts using expected density bins
  mean_scaled_counts_bin1_o1 <- (expected_scaled_counts_o1[1, ] + expected_scaled_counts_o1[2, ]) / 2
  mean_scaled_counts_bin2_o1 <- (expected_scaled_counts_o1[6, ] + expected_scaled_counts_o1[7, ] +
    expected_scaled_counts_o1[8, ]) / 3
  mean_scaled_counts_bin3_o1 <- (expected_scaled_counts_o1[3, ] + expected_scaled_counts_o1[4, ] +
    expected_scaled_counts_o1[5, ] + expected_scaled_counts_o1[9, ] +
    expected_scaled_counts_o1[10, ]) / 5
  expected_scaled_density_binned_counts_o1 <- rbind(
    mean_scaled_counts_bin1_o1, mean_scaled_counts_bin2_o1, mean_scaled_counts_bin3_o1
  )
  rownames(expected_scaled_density_binned_counts_o1) <- 1:3
  # Order 2: Expected interval indexes: 1, 3, 3, 1, 3, 2, 2, 4, 4, 4
  mean_scaled_counts_bin1_o2 <- (expected_scaled_counts_o2[1, ] + expected_scaled_counts_o2[4, ]) / 2
  mean_scaled_counts_bin2_o2 <- (expected_scaled_counts_o2[2, ] + expected_scaled_counts_o2[6, ] +
    expected_scaled_counts_o2[7, ]) / 3
  mean_scaled_counts_bin3_o2 <- (expected_scaled_counts_o2[3, ] + expected_scaled_counts_o2[5, ]) / 2
  mean_scaled_counts_bin4_o2 <- (expected_scaled_counts_o2[8, ] + expected_scaled_counts_o2[9, ] +
    expected_scaled_counts_o2[10, ]) / 3
  expected_scaled_density_binned_counts_o2 <- rbind(
    mean_scaled_counts_bin1_o2, mean_scaled_counts_bin2_o2, mean_scaled_counts_bin3_o2,
    mean_scaled_counts_bin4_o2
  )
  rownames(expected_scaled_density_binned_counts_o2) <- 1:4

  # Calculate scaled binned counts with density_binned_counts (function to test)
  actual_scaled_density_binned_counts_o1 <- density_binned_counts(
    expected_counts_o1,
    expected_interval_indexes_o1,
    agg_fn = agg_fn,
    scale_fn = scale_fn,
    max_graphlet_size = max_graphlet_size
  )

  actual_scaled_density_binned_counts_o2 <- density_binned_counts(
    expected_counts_o2,
    expected_interval_indexes_o2,
    agg_fn = agg_fn,
    scale_fn = scale_fn,
    max_graphlet_size = max_graphlet_size
  )

  # Check actual output vs expected
  expect_equal(
    actual_scaled_density_binned_counts_o1,
    expected_scaled_density_binned_counts_o1
  )
  expect_equal(
    actual_scaled_density_binned_counts_o2,
    expected_scaled_density_binned_counts_o2
  )
})

context("Measures Netdis: Expected graphlet counts")
test_that("netdis_expected_counts_ego works for graphlets up to 4 nodes", {
  # Helper function to generate graphs with known density and number of nodes
  rand_graph <- function(num_nodes, density) {
    max_edges <- choose(num_nodes, 2)
    num_edges <- density * max_edges
    igraph::erdos.renyi.game(num_nodes, num_edges, "gnm",
      loops = FALSE, directed = FALSE
    )
  }
  # Set up some dummy reference density breaks and scaled reference counts
  density_breaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
  scaled_reference_counts <- rbind(
    c(1, 2, 3, 4, 5, 6, 7, 8, 9),
    c(11, 12, 13, 14, 15, 16, 17, 18, 19),
    c(21, 22, 23, 24, 25, 26, 27, 28, 29),
    c(31, 32, 33, 34, 35, 36, 37, 38, 39),
    c(41, 42, 43, 44, 45, 46, 47, 48, 49),
    c(51, 52, 53, 54, 55, 56, 57, 58, 59),
    c(61, 62, 63, 64, 65, 66, 67, 68, 69),
    c(71, 72, 73, 74, 75, 76, 77, 78, 79),
    c(81, 82, 83, 84, 85, 86, 87, 88, 89),
    c(91, 92, 93, 94, 95, 96, 97, 98, 99)
  )
  graphlet_labels <- c("G0", "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8")
  colnames(scaled_reference_counts) <- graphlet_labels
  rownames(scaled_reference_counts) <- 1:10
  graphlet_sizes <- c(2, 3, 3, 4, 4, 4, 4, 4, 4)
  names(graphlet_sizes) <- graphlet_labels
  max_graphlet_size <- 4

  # Generate some test graphs
  densities <- c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95)
  density_indexes <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  num_nodes <- rep(120, 10)
  graphs <- purrr::map2(num_nodes, densities, rand_graph)
  graphlet_counts <- purrr::map(graphs,
    count_graphlets_for_graph,
    max_graphlet_size = max_graphlet_size
  )
  # WITH scale_fn = NULL (bin counts directly with no scaling)
  # Helper function to calculate expected expected graphlet counts
  expected_expected_graphlet_counts_fn <- function(density_index) {
    scaled_reference_counts[density_index, ]
  }
  # Determine expected and actual expected graphlet counts
  expected_expected_graphlet_counts <-
    purrr::map(density_indexes, expected_expected_graphlet_counts_fn)
  actual_expected_graphlet_counts <-
    purrr::map(graphlet_counts, netdis_expected_counts_ego,
      max_graphlet_size = max_graphlet_size,
      density_breaks = density_breaks,
      density_binned_reference_counts = scaled_reference_counts,
      scale_fn = NULL
    )

  # Loop over each graph and compare expected with actual
  # NOTE: v2.0.0 of testthat library made a breaking change that means using
  # map, mapply etc can cause failures under certain conditions
  # See: https://github.com/r-lib/testthat/releases/tag/v2.0.0
  for (i in 1:length(actual_expected_graphlet_counts)) {
    expect_equal(
      actual_expected_graphlet_counts[i],
      expected_expected_graphlet_counts[i]
    )
  }

  # WITH scale_fn = count_graphlet_tuples (default netdis from paper)
  # Helper function to calculate expected expected graphlet counts
  expected_expected_graphlet_counts_fn <- function(density_index, node_count) {
    reference_counts <- scaled_reference_counts[density_index, ]
    reference_counts * choose(node_count, graphlet_sizes)
  }
  # Determine expected and actual expected graphlet counts
  expected_expected_graphlet_counts <-
    purrr::map2(density_indexes, num_nodes, expected_expected_graphlet_counts_fn)
  actual_expected_graphlet_counts <-
    purrr::map(graphlet_counts, netdis_expected_counts_ego,
      max_graphlet_size = max_graphlet_size,
      density_breaks = density_breaks,
      density_binned_reference_counts = scaled_reference_counts,
      scale_fn = count_graphlet_tuples
    )
  # Loop over each graph and compare expected with actual
  # NOTE: v2.0.0 of testthat library made a breaking change that means using
  # map, mapply etc can cause failures under certain conditions
  # See: https://github.com/r-lib/testthat/releases/tag/v2.0.0
  for (i in 1:length(actual_expected_graphlet_counts)) {
    expect_equal(
      actual_expected_graphlet_counts[i],
      expected_expected_graphlet_counts[i]
    )
  }
})

test_that("netdis_expected_counts works for graphlets up to 4 nodes", {
  # Helper function to generate graphs with known density and number of nodes
  # Set up a small sample network with at least one ego-network that contains
  # at least one of each graphlets
  elist <- rbind(
    c("n1", "n2"),
    c("n2", "n3"),
    c("n1", "n4"),
    c("n2", "n5"),
    c("n1", "n6"),
    c("n1", "n7"),
    c("n2", "n4"),
    c("n4", "n6"),
    c("n6", "n8"),
    c("n7", "n8"),
    c("n7", "n9"),
    c("n7", "n10"),
    c("n8", "n9"),
    c("n8", "n10"),
    c("n9", "n10")
  )
  graph <- igraph::graph_from_edgelist(elist, directed = FALSE)
  graphlet_labels <- c("G0", "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8")
  graphlet_sizes <- c(2, 3, 3, 4, 4, 4, 4, 4, 4)
  max_graphlet_size <- 4
  min_ego_edges <- 0
  min_ego_nodes <- 0

  # Get ego network graphlet counts
  graphlet_counts_ego_o1 <- count_graphlets_ego(graph,
    neighbourhood_size = 1,
    min_ego_edges = min_ego_edges,
    min_ego_nodes = min_ego_nodes
  )
  graphlet_counts_ego_o2 <- count_graphlets_ego(graph,
    neighbourhood_size = 2,
    min_ego_edges = min_ego_edges,
    min_ego_nodes = min_ego_nodes
  )
  # Set manually-verified node counts and densities
  # 1. Ego-networks of order 1
  num_nodes_o1 <- c(5, 5, 2, 4, 2, 4, 5, 5, 4, 4)
  num_edges_o1 <- c(6, 5, 1, 5, 1, 4, 7, 7, 6, 6)
  max_edges_o1 <- choose(num_nodes_o1, 2)
  densities_o1 <- num_edges_o1 / max_edges_o1
  # Order 1 densities should be: 0.6000000 0.5000000 1.0000000 0.8333333 1.0000000 0.6666667 0.7000000 0.7000000 1.0000000 1.0000000
  # 2. Ego-networks of order 2
  num_nodes_o2 <- c(10, 7, 5, 8, 5, 8, 8, 7, 6, 6)
  num_edges_o2 <- c(15, 8, 5, 10, 5, 13, 13, 11, 9, 9)
  max_edges_o2 <- choose(num_nodes_o2, 2)
  densities_o2 <- num_edges_o2 / max_edges_o2
  # Order 2 densities should be: 0.3333333 0.3809524 0.5000000 0.3571429 0.5000000 0.4642857 0.4642857 0.5238095 0.6000000 0.6000000
  # Set manually defined density breaks and indexes
  breaks <- c(0, 0.11, 0.21, 0.31, 0.41, 0.51, 0.61, 0.71, 0.81, 0.91, 1.0)
  density_indexes_o1 <- c(6, 5, 10, 9, 10, 7, 7, 7, 10, 10)
  density_indexes_o2 <- c(4, 4, 5, 4, 5, 5, 5, 6, 6, 6)
  # Set dummy reference counts
  scaled_reference_counts <- rbind(
    c(1, 2, 3, 4, 5, 6, 7, 8, 9),
    c(11, 12, 13, 14, 15, 16, 17, 18, 19),
    c(21, 22, 23, 24, 25, 26, 27, 28, 29),
    c(31, 32, 33, 34, 35, 36, 37, 38, 39),
    c(41, 42, 43, 44, 45, 46, 47, 48, 49),
    c(51, 52, 53, 54, 55, 56, 57, 58, 59),
    c(61, 62, 63, 64, 65, 66, 67, 68, 69),
    c(71, 72, 73, 74, 75, 76, 77, 78, 79),
    c(81, 82, 83, 84, 85, 86, 87, 88, 89),
    c(91, 92, 93, 94, 95, 96, 97, 98, 99)
  )
  colnames(scaled_reference_counts) <- graphlet_labels
  expected_dims <- dim(scaled_reference_counts)
  min_ego_nodes <- 3
  min_ego_edges <- 1

  #-------------------------------------------------------
  # With scale_fn = count_graphlet_tuples (default netdis paper)
  #-------------------------------------------------------
  # Helper function to calculate expected expected graphlet counts
  expected_expected_graphlet_counts_fn <- function(density_index, node_count) {
    reference_counts <- scaled_reference_counts[density_index, ]
    reference_counts * choose(node_count, graphlet_sizes)
  }
  # Calculate expected graphlet counts. NOTE: We expect a matrix with graphlet
  # types as columns and ego networks for nodes in graph as rows
  expected_expected_graphlet_counts_ego_o1 <- t(simplify2array(purrr::map2(
    density_indexes_o1, num_nodes_o1, expected_expected_graphlet_counts_fn
  )))
  expected_expected_graphlet_counts_ego_o2 <- t(simplify2array(purrr::map2(
    density_indexes_o2, num_nodes_o2, expected_expected_graphlet_counts_fn
  )))
  # Sanity check for expected output shape. Should be matrix with graphlet types
  # as columns and nodes as rows
  expect_equal(dim(expected_expected_graphlet_counts_ego_o1), expected_dims)
  expect_equal(dim(expected_expected_graphlet_counts_ego_o2), expected_dims)
  # Set column labels to graphlet names
  colnames(expected_expected_graphlet_counts_ego_o1) <- graphlet_labels
  colnames(expected_expected_graphlet_counts_ego_o2) <- graphlet_labels
  # Set row labels to ego network names
  rownames(expected_expected_graphlet_counts_ego_o1) <- rownames(graphlet_counts_ego_o1)
  rownames(expected_expected_graphlet_counts_ego_o2) <- rownames(graphlet_counts_ego_o1)

  # Calculate actual output of function under test
  actual_expected_graphlet_counts_ego_o1 <-
    netdis_expected_counts(
      graphlet_counts_ego_o1,
      breaks,
      scaled_reference_counts,
      max_graphlet_size,
      scale_fn = count_graphlet_tuples
    )
  actual_expected_graphlet_counts_ego_o2 <-
    netdis_expected_counts(
      graphlet_counts_ego_o2,
      breaks,
      scaled_reference_counts,
      max_graphlet_size,
      scale_fn = count_graphlet_tuples
    )

  # Compare actual to expected
  expect_equal(
    actual_expected_graphlet_counts_ego_o1,
    expected_expected_graphlet_counts_ego_o1
  )
  expect_equal(
    actual_expected_graphlet_counts_ego_o2,
    expected_expected_graphlet_counts_ego_o2
  )

  #-------------------------------------------------------
  # With scale_fn = NULL (take reference counts directly)
  #-------------------------------------------------------
  # Helper function to calculate expected expected graphlet counts
  expected_expected_graphlet_counts_fn <- function(density_index) {
    scaled_reference_counts[density_index, ]
  }
  # Calculate expected graphlet counts. NOTE: We expect a matrix with graphlet
  # types as columns and ego networks for nodes in graph as rows
  expected_expected_graphlet_counts_ego_o1 <- t(simplify2array(purrr::map(
    density_indexes_o1, expected_expected_graphlet_counts_fn
  )))
  expected_expected_graphlet_counts_ego_o2 <- t(simplify2array(purrr::map(
    density_indexes_o2, expected_expected_graphlet_counts_fn
  )))
  # Sanity check for expected output shape. Should be matrix with graphlet types
  # as columns and nodes as rows
  expect_equal(dim(expected_expected_graphlet_counts_ego_o1), expected_dims)
  expect_equal(dim(expected_expected_graphlet_counts_ego_o2), expected_dims)
  # Set column labels to graphlet names
  colnames(expected_expected_graphlet_counts_ego_o1) <- graphlet_labels
  colnames(expected_expected_graphlet_counts_ego_o2) <- graphlet_labels
  # Set row labels to ego network names
  rownames(expected_expected_graphlet_counts_ego_o1) <- rownames(graphlet_counts_ego_o1)
  rownames(expected_expected_graphlet_counts_ego_o2) <- rownames(graphlet_counts_ego_o2)

  # Calculate actual output of function under test
  actual_expected_graphlet_counts_ego_o1 <-
    netdis_expected_counts(
      graphlet_counts_ego_o1,
      breaks,
      scaled_reference_counts,
      max_graphlet_size,
      scale_fn = NULL
    )
  actual_expected_graphlet_counts_ego_o2 <-
    netdis_expected_counts(
      graphlet_counts_ego_o2,
      breaks,
      scaled_reference_counts,
      max_graphlet_size,
      scale_fn = NULL
    )

  # Compare actual to expected
  expect_equal(
    actual_expected_graphlet_counts_ego_o1,
    expected_expected_graphlet_counts_ego_o1
  )
  expect_equal(
    actual_expected_graphlet_counts_ego_o2,
    expected_expected_graphlet_counts_ego_o2
  )
})

context("Netdis: Statistic calculation")
test_that("netdis statistic function output matches manually verified result", {
  # arbitrary counts of correct size for graphlets up to size 5
  counts_1 <- c(
    11, 11, 13, 9, 12, 10, 14, 9, 13, 10, 10, 7, 9, 12, 6, 12, 9, 12,
    9, 7, 15, 7, 5, 12, 16, 10, 10, 8, 9, 14
  )
  counts_2 <- c(
    12, 11, 6, 10, 15, 7, 10, 8, 7, 7, 7, 13, 9, 14, 7, 12,
    7, 10, 9, 11, 7, 7, 11, 8, 10, 14, 8, 16, 14, 10
  )

  # add graphlet names
  ids <- graphlet_key(5)$id
  names(counts_1) <- ids
  names(counts_2) <- ids

  # manually verified results
  expected_netdis_3 <- 0.03418796
  expected_netdis_4 <- 0.02091792
  expected_netdis_5 <- 0.03826385

  # check function to test
  actual_netdis_3 <- netdis(counts_1, counts_2, 3)
  actual_netdis_4 <- netdis(counts_1, counts_2, 4)
  actual_netdis_5 <- netdis(counts_1, counts_2, 5)

  expect_equal(expected_netdis_3, actual_netdis_3)
  expect_equal(expected_netdis_4, actual_netdis_4)
  expect_equal(expected_netdis_5, actual_netdis_5)
})
test_that("netdis_uptok gives expected netdis result for graphlets up to size k", {
  # arbitrary counts of correct size for graphlets up to size 5
  counts_1 <- c(
    11, 11, 13, 9, 12, 10, 14, 9, 13, 10, 10, 7, 9, 12, 6, 12, 9, 12,
    9, 7, 15, 7, 5, 12, 16, 10, 10, 8, 9, 14
  )
  counts_2 <- c(
    12, 11, 6, 10, 15, 7, 10, 8, 7, 7, 7, 13, 9, 14, 7, 12,
    7, 10, 9, 11, 7, 7, 11, 8, 10, 14, 8, 16, 14, 10
  )

  # add graphlet names
  ids <- graphlet_key(5)$id
  names(counts_1) <- ids
  names(counts_2) <- ids

  # manually verified results
  expected_netdis <- c(0.03418796, 0.02091792, 0.03826385)
  names(expected_netdis) <- c("netdis3", "netdis4", "netdis5")

  # check function to test
  actual_netdis <- netdis_uptok(counts_1, counts_2, 5)

  expect_equal(expected_netdis, actual_netdis)
})

context("Netdis: full calculation pipeline")
test_that("netdis_many_to_many gives expected result", {
  # Set source directory for Virus PPI graph edge files
  source_dir <- system.file(file.path("extdata", "VRPINS"), package = "netdist")

  # Load query and reference graphs
  graphs <- read_simple_graphs(source_dir, format = "ncol", pattern = "*")
  graphs <- graphs[c("EBV", "ECL", "HSV-1", "KSHV")]

  ref_path <- system.file(file.path("extdata", "random", "ER_1250_10_1"),
    package = "netdist"
  )
  ref_graph <- read_simple_graph(ref_path, format = "ncol")

  # set parameters
  max_graphlet_size <- 4
  neighbourhood_size <- 2
  min_ego_nodes <- 3
  min_ego_edges <- 1

  # manually verified results
  # $netdis
  #              [,1]        [,2]       [,3]      [,4]      [,5]         [,6]
  # netdis3 0.1846655 0.008264222 0.01005385 0.2065762 0.2091241 0.0001335756
  # netdis4 0.1749835 0.165264120 0.01969246 0.2917612 0.2215579 0.0760242643
  #
  # $comp_spec
  # name_a name_b index_a index_b
  # 1    EBV    ECL       1       2
  # 2    EBV  HSV-1       1       3
  # 3    EBV   KSHV       1       4
  # 4    ECL  HSV-1       2       3
  # 5    ECL   KSHV       2       4
  # 6  HSV-1   KSHV       3       4
  expected_netdis_netdis <- matrix(nrow = 2, ncol = 6)
  expected_netdis_netdis[1, ] <- c(
    0.1846655, 0.008264222, 0.01005385,
    0.2065762, 0.2091241, 0.0001335756
  )
  expected_netdis_netdis[2, ] <- c(
    0.1749835, 0.165264120, 0.01969246,
    0.2917612, 0.2215579, 0.0760242643
  )
  rownames(expected_netdis_netdis) <- c("netdis3", "netdis4")

  expected_netdis_comp_spec <- cross_comparison_spec(
    list(
      "EBV" = c(),
      "ECL" = c(),
      "HSV-1" = c(),
      "KSHV" = c()
    )
  )

  expected_netdis <- list(
    netdis = expected_netdis_netdis,
    comp_spec = expected_netdis_comp_spec
  )


  # Calculate netdis statistics
  actual_netdis <- netdis_many_to_many(graphs,
    ref_graph,
    max_graphlet_size = max_graphlet_size,
    neighbourhood_size = neighbourhood_size,
    min_ego_nodes = min_ego_nodes,
    min_ego_edges = min_ego_edges
  )

  # Check results as expected
  expect_equal(expected_netdis, actual_netdis, tolerance = .001, scale = 1)
})

context("Netdis: functions for different pairwise comparisons")
test_that("netdis_one_to_one gives expected result", {
  # Set source directory for Virus PPI graph edge files
  source_dir <- system.file(file.path("extdata", "VRPINS"), package = "netdist")

  # Load query and reference graphs
  graph_1 <- read_simple_graph(file.path(source_dir, "EBV.txt"),
    format = "ncol"
  )

  graph_2 <- read_simple_graph(file.path(source_dir, "ECL.txt"),
    format = "ncol"
  )

  ref_path <- system.file(file.path("extdata", "random", "ER_1250_10_1"),
    package = "netdist"
  )
  ref_graph <- read_simple_graph(ref_path, format = "ncol")

  # set parameters
  max_graphlet_size <- 4
  neighbourhood_size <- 2
  min_ego_nodes <- 3
  min_ego_edges <- 1

  # manually verified results
  expected_netdis <- c(0.1846655, 0.1749835)
  names(expected_netdis) <- c("netdis3", "netdis4")

  # check function to test
  actual_netdis <- netdis_one_to_one(graph_1,
    graph_2,
    ref_graph,
    max_graphlet_size = max_graphlet_size,
    neighbourhood_size = neighbourhood_size,
    min_ego_nodes = min_ego_nodes,
    min_ego_edges = min_ego_edges
  )

  expect_equal(expected_netdis, actual_netdis, tolerance = .001, scale = 1)
})
test_that("netdis_one_to_many gives expected result", {
  # Set source directory for Virus PPI graph edge files
  source_dir <- system.file(file.path("extdata", "VRPINS"), package = "netdist")

  # Load query and reference graphs
  graphs <- read_simple_graphs(source_dir, format = "ncol", pattern = "*")
  graph_1 <- graphs$EBV
  graphs_compare <- graphs[c("ECL", "HSV-1", "KSHV", "VZV")]

  ref_path <- system.file(file.path("extdata", "random", "ER_1250_10_1"),
    package = "netdist"
  )
  ref_graph <- read_simple_graph(ref_path, format = "ncol")

  # set parameters
  max_graphlet_size <- 4
  neighbourhood_size <- 2
  min_ego_nodes <- 3
  min_ego_edges <- 1

  # manually verified results
  #         ECL       HSV-1       KSHV       VZV
  # netdis3 0.1846655 0.008264222 0.01005385 0.006777578
  # netdis4 0.1749835 0.165264120 0.01969246 0.159711160
  expected_netdis <- matrix(nrow = 2, ncol = 4)
  colnames(expected_netdis) <- c("ECL", "HSV-1", "KSHV", "VZV")
  rownames(expected_netdis) <- c("netdis3", "netdis4")
  expected_netdis[1, ] <- c(0.1846655, 0.008264222, 0.01005385, 0.006777578)
  expected_netdis[2, ] <- c(0.1749835, 0.165264120, 0.01969246, 0.159711160)

  # Calculate netdis statistics
  actual_netdis <- netdis_one_to_many(graph_1, graphs_compare,
    ref_graph,
    max_graphlet_size = max_graphlet_size,
    neighbourhood_size = neighbourhood_size,
    min_ego_nodes = min_ego_nodes,
    min_ego_edges = min_ego_edges
  )

  # Check results as expected
  expect_equal(expected_netdis, actual_netdis, tolerance = .001, scale = 1)
})

context("Netdis: error if no query graphs or graphlet counts provided")
test_that("netdis functions error when no query graphs provided", {
  # dummy values to use for other parameters
  ref_graph <- NULL
  max_graphlet_size <- 4
  neighbourhood_size <- 2
  min_ego_nodes <- 3
  min_ego_edges <- 1
  comparisons <- "many-to-many"

  elist <- rbind(
    c("n1", "n2"),
    c("n2", "n3"),
    c("n1", "n4"),
    c("n2", "n5"),
    c("n1", "n6"),
    c("n1", "n7"),
    c("n2", "n4"),
    c("n4", "n6"),
    c("n6", "n8"),
    c("n7", "n8"),
    c("n7", "n9"),
    c("n7", "n10"),
    c("n8", "n9"),
    c("n8", "n10"),
    c("n9", "n10")
  )
  graph <- igraph::graph_from_edgelist(elist, directed = FALSE)

  # one to one function
  expect_error(
    netdis_one_to_one(
      graph_1 = graph,
      graph_2 = NULL,
      ref_graph = ref_graph,
      max_graphlet_size = max_graphlet_size,
      neighbourhood_size = neighbourhood_size,
      min_ego_nodes = min_ego_nodes,
      min_ego_edges = min_ego_edges,
      graphlet_counts_1 = NULL,
      graphlet_counts_2 = NULL
    )
  )
  expect_error(
    netdis_one_to_one(
      graph_1 = NULL,
      graph_2 = graph,
      ref_graph = ref_graph,
      max_graphlet_size = max_graphlet_size,
      neighbourhood_size = neighbourhood_size,
      min_ego_nodes = min_ego_nodes,
      min_ego_edges = min_ego_edges,
      graphlet_counts_1 = NULL,
      graphlet_counts_2 = NULL
    )
  )

  # one to many function
  expect_error(
    netdis_one_to_many(
      graph_1 = graph,
      graphs_compare = NULL,
      ref_graph = ref_graph,
      max_graphlet_size = max_graphlet_size,
      neighbourhood_size = neighbourhood_size,
      min_ego_nodes = min_ego_nodes,
      min_ego_edges = min_ego_edges,
      graphlet_counts_1 = NULL,
      graphlet_counts_compare = NULL
    )
  )
  expect_error(
    netdis_one_to_many(
      graph_1 = NULL,
      graphs_compare = list(graph_1 = graph, graph_2 = graph),
      ref_graph = ref_graph,
      max_graphlet_size = max_graphlet_size,
      neighbourhood_size = neighbourhood_size,
      min_ego_nodes = min_ego_nodes,
      min_ego_edges = min_ego_edges,
      graphlet_counts_1 = NULL,
      graphlet_counts_compare = NULL
    )
  )
  # many to many function
  expect_error(
    netdis_many_to_many(
      graphs = NULL,
      comparisons = comparisons,
      ref_graph = ref_graph,
      max_graphlet_size = max_graphlet_size,
      neighbourhood_size = neighbourhood_size,
      min_ego_nodes = min_ego_nodes,
      min_ego_edges = min_ego_edges,
      graphlet_counts = NULL
    )
  )
})

context("Netdis: constant expected counts")
test_that("netdis_centred_graphlet_counts correctly interprets numeric
          ref_binned_graphlet_counts value", {
  # dummy counts
  graphlet_counts <- rbind(
    c(1, 2, 3, 4, 5, 6, 7, 8, 9),
    c(11, 12, 13, 14, 15, 16, 17, 18, 19),
    c(21, 22, 23, 24, 25, 26, 27, 28, 29),
    c(31, 32, 33, 34, 35, 36, 37, 38, 39),
    c(41, 42, 43, 44, 45, 46, 47, 48, 49),
    c(51, 52, 53, 54, 55, 56, 57, 58, 59),
    c(61, 62, 63, 64, 65, 66, 67, 68, 69),
    c(71, 72, 73, 74, 75, 76, 77, 78, 79),
    c(81, 82, 83, 84, 85, 86, 87, 88, 89),
    c(91, 92, 93, 94, 95, 96, 97, 98, 99)
  )
  graphlet_labels <- c("G0", "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8")
  colnames(graphlet_counts) <- graphlet_labels
  rownames(graphlet_counts) <- 1:10
  max_graphlet_size <- 4

  # netdis_centred_graphlet_counts with ref_binned_graphlet_counts=0 should
  # perform no centring, i.e. centred_counts should equal input graphlet_counts
  centred_counts <- netdis_centred_graphlet_counts(
    graphlet_counts = graphlet_counts,
    ref_ego_density_bins = NULL,
    ref_binned_graphlet_counts = 0,
    binning_fn = NULL,
    bin_counts_fn = NULL,
    exp_counts_fn = NULL,
    max_graphlet_size = 4
  )

  expect_equal(centred_counts, graphlet_counts)
})

context("Netdis: Geometric Poisson Approximation")
test_that("netdis_one_to_one gives expected result when using geometric Poisson
          approximation", {
  # Set source directory for Virus PPI graph edge files
  source_dir <- system.file(file.path("extdata", "VRPINS"), package = "netdist")

  # Load query and reference graphs
  graph_1 <- read_simple_graph(file.path(source_dir, "EBV.txt"),
    format = "ncol"
  )

  graph_2 <- read_simple_graph(file.path(source_dir, "ECL.txt"),
    format = "ncol"
  )

  # set parameters
  max_graphlet_size <- 4
  neighbourhood_size <- 2
  min_ego_nodes <- 3
  min_ego_edges <- 1

  # manually verified result for graphlets of size 4
  # verified using a different implementation of geometric poisson with these
  # networks.
  expected_netdis4 <- 0.1892716

  # check function to test
  actual_netdis <- netdis_one_to_one(graph_1,
    graph_2,
    ref_graph = NULL,
    max_graphlet_size = max_graphlet_size,
    neighbourhood_size = neighbourhood_size,
    min_ego_nodes = min_ego_nodes,
    min_ego_edges = min_ego_edges,
  )

  expect_equal(expected_netdis4, actual_netdis[["netdis4"]],
    tolerance = .0001, scale = 1
  )
})

context("Netdis: error if max_graphlet_size is not 3, 4 or 5")
test_that("netdis_uptok errors for unsupported max_graphlet_size", {
  # dummy counts values
  counts_1 <- c(
    11, 11, 13, 9, 12, 10, 14, 9, 13, 10, 10, 7, 9, 12, 6, 12, 9, 12,
    9, 7, 15, 7, 5, 12, 16, 10, 10, 8, 9, 14
  )
  counts_2 <- c(
    12, 11, 6, 10, 15, 7, 10, 8, 7, 7, 7, 13, 9, 14, 7, 12,
    7, 10, 9, 11, 7, 7, 11, 8, 10, 14, 8, 16, 14, 10
  )
  ids <- graphlet_key(5)$id
  names(counts_1) <- ids
  names(counts_2) <- ids

  # graphlet size greater than 5
  expect_error(netdis_uptok(counts_1, counts_2, 6))

  #  graphlet size less than 3
  expect_error(netdis_uptok(counts_1, counts_2, 2))
})

context("Netdis: works correctly when using a single density bin")
test_that("netdis single density bin works correctly", {
  # TODO
})
