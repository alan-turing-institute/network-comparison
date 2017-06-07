
context("Measures Netdis: Adaptive binning")
test_that("adaptive_breaks merges 2 lowest bins where only first bin is below minimum", {
  min_count <- 5
  x <- c(1.5, rep(2.2, min_count), rep(3.5, min_count), rep(4.5, min_count), 
         rep(5.5, min_count), rep(6.5, min_count + 1))
  initial_breaks <- 1:7
  final_breaks_actual <- adaptive_breaks(x, min_count = min_count, breaks = initial_breaks)
  final_breaks_expected <- c(1, 3, 4, 5, 6, 7)
  
  expect_equal(final_breaks_actual, final_breaks_expected)
})

test_that("adaptive_breaks merges 3 lowest bins where lowest 2 combined are below minimum", {
  min_count <- 5
  x <- c(1.5, rep(2.2, 2), rep(3.5, min_count), rep(4.5, min_count), 
         rep(5.5, min_count), rep(6.5, min_count + 1))
  initial_breaks <- 1:7
  final_breaks_actual <- adaptive_breaks(x, min_count = min_count, breaks = initial_breaks)
  final_breaks_expected <- c(1, 4, 5, 6, 7)
  
  expect_equal(final_breaks_actual, final_breaks_expected)
})

test_that("adaptive_breaks merges pair of bins in middle", {
  min_count <- 5
  x <- c(rep(1.6, min_count), rep(2.2, min_count), rep(3.5, 2), rep(4.5, 3), 
         rep(5.5, min_count), rep(6.5, min_count + 1))
  initial_breaks <- 1:7
  final_breaks_actual <- adaptive_breaks(x, min_count = min_count, breaks = initial_breaks)
  final_breaks_expected <- c(1, 2, 3, 5, 6, 7)
  
  expect_equal(final_breaks_actual, final_breaks_expected)
})

test_that("adaptive_breaks merges two spearated pairs of bins in middle", {
  min_count <- 5
  x <- c(rep(1.6, min_count), rep(2.2, 2), rep(3.5, 3), rep(4.5, min_count), 
         rep(5.5, 3), rep(6.5, 2), rep(7.8, min_count))
  initial_breaks <- 1:8
  final_breaks_actual <- adaptive_breaks(x, min_count = min_count, breaks = initial_breaks)
  final_breaks_expected <- c(1, 2, 4, 5, 7, 8)
  
  expect_equal(final_breaks_actual, final_breaks_expected)
})

test_that("adaptive_breaks merges 2 uppermost bins where both are below minimum", {
  min_count <- 5
  x <- c(rep(1.5, min_count), rep(2.2, min_count), rep(3.5, min_count),
         rep(4.5, min_count), rep(5.5, 2), rep(6.5, 3))
  initial_breaks <- 1:7
  final_breaks_actual <- adaptive_breaks(x, min_count = min_count, breaks = initial_breaks)
  final_breaks_expected <- c(1, 2,3, 4, 5, 7)

  expect_equal(final_breaks_actual, final_breaks_expected)
})

test_that("adaptive_breaks merges 2 uppermost bins where only last bin is below minimum", {
  min_count <- 5
  x <- c(rep(1.5, min_count), rep(2.2, min_count), rep(3.5, min_count),
         rep(4.5, min_count), rep(5.5, min_count), rep(6.5, 3))
  initial_breaks <- 1:7
  final_breaks_actual <- adaptive_breaks(x, min_count = min_count, breaks = initial_breaks)
  final_breaks_expected <- c(1, 2, 3, 4, 5, 7)
  
  expect_equal(final_breaks_actual, final_breaks_expected)
})

test_that("adaptive_breaks merges bins with no members with the next bin", {
  min_count <- 5
  x <- c(rep(1.5, min_count), rep(5.5, min_count), rep(6.5, min_count))
  initial_breaks <- 1:7
  final_breaks_actual <- adaptive_breaks(x, min_count = min_count, breaks = initial_breaks)
  final_breaks_expected <- c(1, 2, 6, 7)
  
  expect_equal(final_breaks_actual, final_breaks_expected)
})


test_that("adaptive_breaks merges 2 bins below minimum, plus the empty bins between them", {
  min_count <- 5
  x <- c(rep(1.5, min_count), rep(2.3, 1), rep(5.5, 4), rep(6.5, min_count))
  initial_breaks <- 1:7
  final_breaks_actual <- adaptive_breaks(x, min_count = min_count, breaks = initial_breaks)
  final_breaks_expected <- c(1, 2, 6, 7)
  
  expect_equal(final_breaks_actual, final_breaks_expected)
})

context("Measures Netdis:  Adaptively binned densities")
test_that("binned_densities_adaptive works", {
  # Helper function
  test_binning <- function(densities, min_counts_per_interval, num_intervals,
                           expected_breaks, expected_interval_indexes) {
    # Set up expected output
    expected <- list(densities = densities, 
                     interval_indexes = expected_interval_indexes, 
                     breaks = expected_breaks)
    # Calculate actual output
    actual <- binned_densities_adaptive(
      densities, min_counts_per_interval = min_counts_per_interval,
      num_intervals = num_intervals)
    # Check actual matches expected
    expect_equal(actual, expected)
  }
  # Test 1:
  densities <- c(0, 0.099, 0.2, 0.299, 0.4, 0.49, 0.6, 0.699, 0.8, 0.899, 1.0)
  min_counts_per_interval <- 2
  num_intervals <- 100
  expected_breaks <-c(0, 0.1, 0.3, 0.5, 0.7, 1.0)
  expected_interval_indexes <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5)
  test_binning(densities, min_counts_per_interval = min_counts_per_interval,
               num_intervals = num_intervals, expected_breaks = expected_breaks,
               expected_interval_indexes = expected_interval_indexes)
  # Test 2:
  densities <- c(0, 0.012, 0.099, 0.201, 0.299, 0.402, 0.49, 0.596, 0.699, 0.803, 0.899, 1.0)
  min_counts_per_interval <- 2
  num_intervals <- 100
  expected_breaks <-c(0, 0.02, 0.21, 0.41, 0.6, 0.81, 1.0)
  expected_interval_indexes <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6)
  test_binning(densities, min_counts_per_interval = min_counts_per_interval,
               num_intervals = num_intervals, expected_breaks = expected_breaks,
               expected_interval_indexes = expected_interval_indexes)
})

context("Measures Netdis: Graphlet tuples")
test_message <-
  paste("count_graphlet_tuples and count_graphlet_tuples_ego give",
        "choose(node_count, graphlet_size) for each graph + graphlet",
        "combination",sep = " ")
test_that(test_message, {
  # Create some test graphs with known node counts (this is the only graph 
  # property we care about for this test)
  graph_n11 <- igraph::erdos.renyi.game(11, p = 1, type = "gnp")
  graph_n37 <- igraph::erdos.renyi.game(37, p = 1, type = "gnp")
  graph_n73 <- igraph::erdos.renyi.game(73, p = 1, type = "gnp")
  # Calculate expected graph tuple count for graphlets of various sizes. There
  # is 1 graphlet of size 1, 2 of size 3, 6 of size 4, and 21 of size 5
  graphlet_tuple_counts <- function(n, max_graphlet_size) {
    if(max_graphlet_size >= 2) {
      tuple_counts <- rep(choose(n, 2), 1)
    }
    if(max_graphlet_size >= 3) {
      tuple_counts <- c(tuple_counts, rep(choose(n, 3), 2))
    }
    if(max_graphlet_size >= 4) {
      tuple_counts <- c(tuple_counts, rep(choose(n, 4), 6))
    }
    if(max_graphlet_size >= 5) {
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
  actual_tuple_count_n11_gs4 <- count_graphlet_tuples(graph_n11, 4)
  actual_tuple_count_n37_gs4 <- count_graphlet_tuples(graph_n37, 4)
  actual_tuple_count_n73_gs4 <- count_graphlet_tuples(graph_n73, 4)
  actual_tuple_count_n11_gs5 <- count_graphlet_tuples(graph_n11, 5)
  actual_tuple_count_n37_gs5 <- count_graphlet_tuples(graph_n37, 5)
  actual_tuple_count_n73_gs5 <- count_graphlet_tuples(graph_n73, 5)
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
  graphlet_tuple_counts_ego <- function(ego_networks, max_graphlet_size) {
    t(sapply(ego_networks, FUN = function(g) {
      graphlet_tuple_counts(length(igraph::V(g)), max_graphlet_size)}))
  }
  # Generate ego networks for each graph
  graph_n11_ego1 <- make_named_ego_graph(graph_n11, order = 1)
  graph_n37_ego1 <- make_named_ego_graph(graph_n37, order = 1)
  graph_n73_ego1 <- make_named_ego_graph(graph_n73, order = 1)
  graph_n11_ego2 <- make_named_ego_graph(graph_n11, order = 2)
  graph_n37_ego2 <- make_named_ego_graph(graph_n37, order = 2)
  graph_n73_ego2 <- make_named_ego_graph(graph_n73, order = 2)
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
  # 1. For ego-networks of order 1
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

context("Measures Netdis: Ego-network scaled graphlet outputs for manually verified networks")
test_that("Ego-network 4-node graphlet counts match manually verified totals",{
  # Set up a small sample network with at least one ego-network that contains
  # at least one of each graphlets
  elist <- rbind(
    c("n1","n2"),
    c("n2","n3"),
    c("n1","n4"),
    c("n2","n5"),
    c("n1","n6"),
    c("n1","n7"),
    c("n2","n4"),
    c("n4","n6"),
    c("n6","n8"),
    c("n7","n8"),
    c("n7","n9"),
    c("n7","n10"),
    c("n8","n9"),
    c("n8","n10"),
    c("n9","n10")
  )
  graph <- igraph::graph_from_edgelist(elist, directed = FALSE)
  
  # Set node and graphlet labels to use for row and col names in expected counts
  node_labels <- igraph::V(graph)$name
  graphlet_labels <- c("G0", "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8")
  
  # Count graphlets in each ego network of the graph with neighbourhood sizes of 1 and 2
  max_graphlet_size <- 4
  actual_counts_order_1 <- 
    count_graphlets_ego_scaled(graph, max_graphlet_size = max_graphlet_size, 
                        neighbourhood_size = 1)
  actual_counts_order_2 <- 
    count_graphlets_ego_scaled(graph, max_graphlet_size = max_graphlet_size, 
                        neighbourhood_size = 2)
  graphlet_key <- graphlet_key(max_graphlet_size)
  k <- graphlet_key$node_count
  # Set manually verified counts
  # 1-step ego networks
  expected_counts_order_1 <- rbind(
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
  rownames(expected_counts_order_1) <- node_labels
  colnames(expected_counts_order_1) <- graphlet_labels
  # 2-step ego networks
  expected_counts_order_2 <- rbind(
    c(15, 18, 6, 21, 3, 1, 11, 1, 1) / zeros_to_ones(choose(10 , k)),
    c( 8, 10, 2,  6, 3, 0,  4, 1, 0) / zeros_to_ones(choose(7 , k)),
    c( 5,  5, 1,  0, 2, 0,  2, 0, 0) / zeros_to_ones(choose(5 , k)),
    c(10, 14, 2, 11, 3, 1,  5, 1, 0) / zeros_to_ones(choose(8 , k)),
    c( 5,  5, 1,  0, 2, 0,  2, 0, 0) / zeros_to_ones(choose(5 , k)),
    c(13, 13, 6, 15, 1, 1,  9, 1, 1) / zeros_to_ones(choose(8 , k)),
    c(13, 13, 6, 15, 1, 1,  9, 1, 1) / zeros_to_ones(choose(8 , k)),
    c(11, 10, 5, 10 ,0 ,1,  8, 0, 1) / zeros_to_ones(choose(7 , k)),
    c( 9,  8, 4,  4, 0, 1,  6, 0, 1) / zeros_to_ones(choose(6 , k)),
    c( 9,  8, 4,  4, 0, 1,  6, 0, 1) / zeros_to_ones(choose(6 , k))
  )
  rownames(expected_counts_order_2) <- node_labels
  colnames(expected_counts_order_2) <- graphlet_labels
  
  # Test that actual counts match expected with only counts requested (default)
  expect_equal(actual_counts_order_1, expected_counts_order_1)
  expect_equal(actual_counts_order_2, expected_counts_order_2)
  
  # Test that actual counts and returned ego networks match expected
  # 1. Define expected
  expected_ego_networks_order_1 <- make_named_ego_graph(graph, order = 1)
  expected_ego_networks_order_2 <- make_named_ego_graph(graph, order = 2)
  expected_counts_with_networks_order_1 <-
    list(graphlet_counts = expected_counts_order_1,
         ego_networks = expected_ego_networks_order_1)
  expected_counts_with_networks_order_2 <- 
    list(graphlet_counts = expected_counts_order_2,
         ego_networks = expected_ego_networks_order_2)
  # 2. Calculate actual
  actual_counts_with_networks_order_1 <- 
    count_graphlets_ego_scaled(graph, max_graphlet_size = max_graphlet_size, 
                        neighbourhood_size = 1, return_ego_networks = TRUE)
  actual_counts_with_networks_order_2 <- 
    count_graphlets_ego_scaled(graph, max_graphlet_size = max_graphlet_size, 
                        neighbourhood_size = 2, return_ego_networks = TRUE)
  
  # 3. Compare
  # Comparison is not implemented for igraph objects, so convert all igraphs to 
  # indexed edge list and then compare. Do in-situ replacement of igraphs with
  # indexed edge lists to ensure we are checking full properties of returned
  # objects (i.e. named lists with matching elements).
  # 3a. Convert expected and actual ego networks from igraphs to indexed edges
  expected_counts_with_networks_order_1$ego_networks <- 
    purrr::map(expected_counts_with_networks_order_1$ego_networks, 
               graph_to_indexed_edges)
  expected_counts_with_networks_order_2$ego_networks <- 
    purrr::map(expected_counts_with_networks_order_2$ego_networks, 
               graph_to_indexed_edges)
  actual_counts_with_networks_order_1$ego_networks <- 
    purrr::map(actual_counts_with_networks_order_1$ego_networks, 
               graph_to_indexed_edges)
  actual_counts_with_networks_order_2$ego_networks <- 
    purrr::map(actual_counts_with_networks_order_2$ego_networks, 
               graph_to_indexed_edges)
  # 3b. Do comparison
  expect_equal(actual_counts_with_networks_order_1, 
               expected_counts_with_networks_order_1)
  expect_equal(actual_counts_with_networks_order_2, 
               expected_counts_with_networks_order_2)
})

context("Measures Netdis: Ego-network density-binned reference counts for manually verified networks")
test_that("Ego-network 4-node density-binned reference counts match manually verified totals",{
  # Set up a small sample network with at least one ego-network that contains
  # at least one of each graphlets
  elist <- rbind(
    c("n1","n2"),
    c("n2","n3"),
    c("n1","n4"),
    c("n2","n5"),
    c("n1","n6"),
    c("n1","n7"),
    c("n2","n4"),
    c("n4","n6"),
    c("n6","n8"),
    c("n7","n8"),
    c("n7","n9"),
    c("n7","n10"),
    c("n8","n9"),
    c("n8","n10"),
    c("n9","n10")
  )
  graph <- igraph::graph_from_edgelist(elist, directed = FALSE)
  
  # Set parameters for test
  max_graphlet_size = 4
  min_counts_per_interval <- 2
  num_intervals <- 100
  
  # Set node and graphlet labels to use for row and col names in expected counts
  node_labels <- igraph::V(graph)$name
  graphlet_labels <- c("G0", "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8")

  # Set manually verified ego-network node counts and edge densities
  #1 . Ego-networks of order 1
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
    expected_densities_o1, min_counts_per_interval = min_counts_per_interval,
    num_intervals = num_intervals)
  expect_equal(actual_binned_densities_o1, expected_binned_densities_o1)
  # 2. Ego-networks of order 2
  expected_min_break_o2 <- 1/3
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
    expected_densities_o2, min_counts_per_interval = min_counts_per_interval,
    num_intervals = num_intervals)
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
    c(15, 18, 6, 21, 3, 1, 11, 1, 1) / zeros_to_ones(choose(10 , k)),
    c( 8, 10, 2,  6, 3, 0,  4, 1, 0) / zeros_to_ones(choose(7 , k)),
    c( 5,  5, 1,  0, 2, 0,  2, 0, 0) / zeros_to_ones(choose(5 , k)),
    c(10, 14, 2, 11, 3, 1,  5, 1, 0) / zeros_to_ones(choose(8 , k)),
    c( 5,  5, 1,  0, 2, 0,  2, 0, 0) / zeros_to_ones(choose(5 , k)),
    c(13, 13, 6, 15, 1, 1,  9, 1, 1) / zeros_to_ones(choose(8 , k)),
    c(13, 13, 6, 15, 1, 1,  9, 1, 1) / zeros_to_ones(choose(8 , k)),
    c(11, 10, 5, 10 ,0 ,1,  8, 0, 1) / zeros_to_ones(choose(7 , k)),
    c( 9,  8, 4,  4, 0, 1,  6, 0, 1) / zeros_to_ones(choose(6 , k)),
    c( 9,  8, 4,  4, 0, 1,  6, 0, 1) / zeros_to_ones(choose(6 , k))
  )
  rownames(expected_counts_o2) <- node_labels
  colnames(expected_counts_o2) <- graphlet_labels
  
  # Calculate binned average expected counts based on manually verified counts
  # and density bins
  # Order 1: Expected interval indexes: 1, 1, 3, 3, 3, 2, 2, 2, 3, 3
  mean_counts_bin1_o1 <- (expected_counts_o1[1,] + expected_counts_o1[2,]) / 2
  mean_counts_bin2_o1 <- (expected_counts_o1[6,] + expected_counts_o1[7,] + 
                            expected_counts_o1[8,]) / 3
  mean_counts_bin3_o1 <- (expected_counts_o1[3,] + expected_counts_o1[4,] + 
                            expected_counts_o1[5,] + expected_counts_o1[9,] +
                            expected_counts_o1[10,]) / 5
  expected_mean_density_binned_counts_o1 <- rbind(
    mean_counts_bin1_o1, mean_counts_bin2_o1, mean_counts_bin3_o1
  )
  rownames(expected_mean_density_binned_counts_o1) <- 1:3
  # Order 2: Expected interval indexes: 1, 3, 3, 1, 3, 2, 2, 4, 4, 4
  mean_counts_bin1_o2 <- (expected_counts_o2[1,] + expected_counts_o2[4,]) / 2
  mean_counts_bin2_o2 <- (expected_counts_o2[2,] + expected_counts_o2[6,] + 
                            expected_counts_o2[7,]) / 3
  mean_counts_bin3_o2 <- (expected_counts_o2[3,] + expected_counts_o2[5,]) / 2
  mean_counts_bin4_o2 <- (expected_counts_o2[8,] + expected_counts_o2[9,] + 
                            expected_counts_o2[10,]) / 3
  expected_mean_density_binned_counts_o2 <- rbind(
    mean_counts_bin1_o2, mean_counts_bin2_o2, mean_counts_bin3_o2, 
    mean_counts_bin4_o2
  )
  rownames(expected_mean_density_binned_counts_o2) <- 1:4
  
  # Calculate actual output of function under test
  actual_mean_density_binned_counts_o1 <- mean_density_binned_graphlet_counts(
      expected_counts_o1, expected_interval_indexes_o1)
  actual_mean_density_binned_counts_o2 <- mean_density_binned_graphlet_counts(
    expected_counts_o2, expected_interval_indexes_o2)
  
  # Check actual output vs expected
  expect_equal(actual_mean_density_binned_counts_o1, 
               expected_mean_density_binned_counts_o1)
  expect_equal(actual_mean_density_binned_counts_o2, 
               expected_mean_density_binned_counts_o2)
})

context("Measures Netdis: Expected graphlet counts")
test_that("netdis_expected_graphlet_counts works for graphlets up to 4 nodes", {
  # Helper function to generate graphs with known density and number of nodes
  rand_graph <- function(num_nodes, density) {
    max_edges <- choose(num_nodes, 2)
    num_edges <- density * max_edges
    igraph::erdos.renyi.game(num_nodes, num_edges , "gnm", 
                             loops = FALSE, directed = FALSE)
  }
  # Set up some dummy reference density breaks and scaled reference counts
  density_breaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
  scaled_reference_counts <- rbind(
    c( 1,  2,  3,  4,  5,  6,  7,  8,  9),
    c(11, 12, 13, 14, 15, 16, 17, 18, 19),
    c(21, 22, 23, 24, 25, 26, 27, 28, 29),
    c(31, 32, 33, 34, 35, 36, 37, 38, 39),
    c(41, 42, 43, 44, 45, 46, 47, 48, 49),
    c(51, 52, 53, 54, 55, 56, 57, 58, 59),
    c(61, 62, 63, 64, 65, 66, 67, 68 ,69),
    c(71, 72, 73, 74, 75, 76, 77, 78, 79),
    c(81, 82, 83, 84, 85, 86 ,87, 88, 89),
    c(91, 92, 93, 94, 95, 96, 97, 98, 99)
  )
  graphlet_labels <- c("G0", "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8")
  colnames(scaled_reference_counts) <- graphlet_labels
  rownames(scaled_reference_counts) <- 1:10
  graphlet_sizes <- c(2, 3, 3, 4, 4, 4, 4, 4, 4)
  names(graphlet_sizes) <- graphlet_labels
  max_graphlet_size = 4
  
  # Generate some test graphs
  densities <- c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95)
  density_indexes <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  num_nodes <- rep(120, 10)
  graphs <- purrr::map2(num_nodes, densities, rand_graph)
  
  # Helper function to calculate expected expected graphlet counts
  expected_expected_graphlet_counts_fn <- function(density_index, node_count) {
    reference_counts <- scaled_reference_counts[density_index,]
    reference_counts * choose(node_count, graphlet_sizes)
    
  }
  # Determine expected and actual expected graphlet counts
  expected_expected_graphlet_counts <- 
    purrr::map2(density_indexes, num_nodes, expected_expected_graphlet_counts_fn)
  actual_expected_graphlet_counts <- 
    purrr::map(graphs, netdis_expected_graphlet_counts, 
               max_graphlet_size = max_graphlet_size, 
               density_breaks = density_breaks, 
               density_binned_reference_counts = scaled_reference_counts)
  # Map over each graph and compare expected with actual
  purrr::map2(actual_expected_graphlet_counts,
              expected_expected_graphlet_counts, expect_equal)
})

test_that("netdis_expected_graphlet_counts_ego works for graphlets up to 4 nodes", {
  # Helper function to generate graphs with known density and number of nodes
  # Set up a small sample network with at least one ego-network that contains
  # at least one of each graphlets
  elist <- rbind(
    c("n1","n2"),
    c("n2","n3"),
    c("n1","n4"),
    c("n2","n5"),
    c("n1","n6"),
    c("n1","n7"),
    c("n2","n4"),
    c("n4","n6"),
    c("n6","n8"),
    c("n7","n8"),
    c("n7","n9"),
    c("n7","n10"),
    c("n8","n9"),
    c("n8","n10"),
    c("n9","n10")
  )
  graph <- igraph::graph_from_edgelist(elist, directed = FALSE)
  graphlet_labels <- c("G0", "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8")
  graphlet_sizes <- c(2, 3, 3, 4, 4, 4, 4, 4, 4)
  max_graphlet_size = 4
  # Make graph ego networks
  ego_networks_o1 <- make_named_ego_graph(graph, order = 1)
  ego_networks_o2 <- make_named_ego_graph(graph, order = 2)
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
    c( 1,  2,  3,  4,  5,  6,  7,  8,  9),
    c(11, 12, 13, 14, 15, 16, 17, 18, 19),
    c(21, 22, 23, 24, 25, 26, 27, 28, 29),
    c(31, 32, 33, 34, 35, 36, 37, 38, 39),
    c(41, 42, 43, 44, 45, 46, 47, 48, 49),
    c(51, 52, 53, 54, 55, 56, 57, 58, 59),
    c(61, 62, 63, 64, 65, 66, 67, 68 ,69),
    c(71, 72, 73, 74, 75, 76, 77, 78, 79),
    c(81, 82, 83, 84, 85, 86 ,87, 88, 89),
    c(91, 92, 93, 94, 95, 96, 97, 98, 99)
  )
  expected_dims <- dim(scaled_reference_counts)
  min_ego_nodes = 3
  min_ego_edges = 1
  
  # Helper function to calculate expected expected graphlet counts
  expected_expected_graphlet_counts_fn <- function(density_index, node_count) {
    reference_counts <- scaled_reference_counts[density_index,]
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
  rownames(expected_expected_graphlet_counts_ego_o1) <- names(ego_networks_o1)
  rownames(expected_expected_graphlet_counts_ego_o2) <- names(ego_networks_o2)
  # Drop rows for nodes with ewer than minumum required nodes and edges in ego 
  # network
  expected_expected_graphlet_counts_ego_o1 <-
    expected_expected_graphlet_counts_ego_o1[
      (num_nodes_o1 >= min_ego_nodes) & (num_edges_o1 >= min_ego_edges),
    ]
  expected_expected_graphlet_counts_ego_o2 <-
    expected_expected_graphlet_counts_ego_o2[
      (num_nodes_o2 >= min_ego_nodes) & (num_edges_o2 >= min_ego_edges),
      ]
  
  # Calculate actual output of function under test
  actual_expected_graphlet_counts_ego_o1 <- 
    netdis_expected_graphlet_counts_ego(
      graph, max_graphlet_size = max_graphlet_size, 
      neighbourhood_size = 1, density_breaks = breaks, 
      density_binned_reference_counts = scaled_reference_counts,
      min_ego_nodes = min_ego_nodes, min_ego_edges = min_ego_edges)
  actual_expected_graphlet_counts_ego_o2 <- 
    netdis_expected_graphlet_counts_ego(
      graph, max_graphlet_size = max_graphlet_size, 
      neighbourhood_size = 2, density_breaks = breaks, 
      density_binned_reference_counts = scaled_reference_counts,
      min_ego_nodes = min_ego_nodes, min_ego_edges = min_ego_edges)
  
  # Compare actual to expected
  expect_equal(actual_expected_graphlet_counts_ego_o1, 
               actual_expected_graphlet_counts_ego_o1)
  expect_equal(actual_expected_graphlet_counts_ego_o2, 
               expected_expected_graphlet_counts_ego_o2)
})

test_that("netdis_expected_graphlet_counts_ego_fn works for graphlets up to 4 nodes", {
  # Set up a small sample network with at least one ego-network that contains
  # at least one of each graphlets
  elist <- rbind(
    c("n1","n2"),
    c("n2","n3"),
    c("n1","n4"),
    c("n2","n5"),
    c("n1","n6"),
    c("n1","n7"),
    c("n2","n4"),
    c("n4","n6"),
    c("n6","n8"),
    c("n7","n8"),
    c("n7","n9"),
    c("n7","n10"),
    c("n8","n9"),
    c("n8","n10"),
    c("n9","n10")
  )
  graph <- igraph::graph_from_edgelist(elist, directed = FALSE)
  graphlet_labels <- c("G0", "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8")
  graphlet_sizes <- c(2, 3, 3, 4, 4, 4, 4, 4, 4)
  names(graphlet_sizes) <- graphlet_labels
  max_graphlet_size = 4
  # Make graph ego networks
  ego_networks_o1 <- make_named_ego_graph(graph, order = 1)
  ego_networks_o2 <- make_named_ego_graph(graph, order = 2)
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
  # Set manually determined density breaks and indexes, based on a min bin count 
  # of 2 and an initial request for 100 bins
  min_bin_count = 2
  num_bins = 100
  num_breaks = num_bins + 1
  min_density_o1 <- 0.5
  max_density_o1 <- 1.0
  breaks_o1 <-  seq(min_density_o1, max_density_o1,length.out = num_breaks)[c(1, 22, 42, 101)]
  density_indexes_o1 <- c(1, 1, 3, 3, 3, 2, 2, 2, 3, 3)
  min_density_o2 <- 1/3
  max_density_o2 <- 0.6
  breaks_o2 <- seq(min_density_o2, max_density_o2,length.out = num_breaks)[c(1, 10, 51, 64, 101)]
  density_indexes_o2 <- c(1, 2, 3, 1, 3, 2, 2, 4, 4, 4)
  # Guard against errors in manually determined breaks and indexes by checking 
  # against already tested code. This also lets us ensure we handle densities
  # falling exactly on a bin boundary the same as the code under test.
  comp_binned_densities_o1 <- binned_densities_adaptive(
    densities_o1, min_counts_per_interval = min_bin_count, 
    num_intervals = num_bins)
  comp_binned_densities_o2 <- binned_densities_adaptive(
    densities_o2, min_counts_per_interval = min_bin_count, 
    num_intervals = num_bins)
  expect_equal(comp_binned_densities_o1, 
               list(densities = densities_o1,
                    interval_indexes = density_indexes_o1,
                    breaks = breaks_o1))
  expect_equal(comp_binned_densities_o2, 
               list(densities = densities_o2,
                    interval_indexes = density_indexes_o2,
                    breaks = breaks_o2))
  
  # Set manually verified scaled ego-network graphlet counts
  graphlet_key <- graphlet_key(max_graphlet_size)
  k <- graphlet_key$node_count
  # 1-step ego networks
  scaled_reference_counts_o1 <- rbind(
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
  # 2-step ego networks
  scaled_reference_counts_o2 <- rbind(
    c(15, 18, 6, 21, 3, 1, 11, 1, 1) / zeros_to_ones(choose(10 , k)),
    c( 8, 10, 2,  6, 3, 0,  4, 1, 0) / zeros_to_ones(choose(7 , k)),
    c( 5,  5, 1,  0, 2, 0,  2, 0, 0) / zeros_to_ones(choose(5 , k)),
    c(10, 14, 2, 11, 3, 1,  5, 1, 0) / zeros_to_ones(choose(8 , k)),
    c( 5,  5, 1,  0, 2, 0,  2, 0, 0) / zeros_to_ones(choose(5 , k)),
    c(13, 13, 6, 15, 1, 1,  9, 1, 1) / zeros_to_ones(choose(8 , k)),
    c(13, 13, 6, 15, 1, 1,  9, 1, 1) / zeros_to_ones(choose(8 , k)),
    c(11, 10, 5, 10 ,0 ,1,  8, 0, 1) / zeros_to_ones(choose(7 , k)),
    c( 9,  8, 4,  4, 0, 1,  6, 0, 1) / zeros_to_ones(choose(6 , k)),
    c( 9,  8, 4,  4, 0, 1,  6, 0, 1) / zeros_to_ones(choose(6 , k))
  )
  min_ego_nodes = 3
  min_ego_edges = 1
  # Drop rows for nodes with ewer than minumum required nodes and edges in ego 
  # network
  scaled_reference_counts_o1 <-
    scaled_reference_counts_o1[
      (num_nodes_o1 >= min_ego_nodes) & (num_edges_o1 >= min_ego_edges),
      ]
  scaled_reference_counts_o2 <-
    scaled_reference_counts_o2[
      (num_nodes_o2 >= min_ego_nodes) & (num_edges_o2 >= min_ego_edges),
      ]
  density_indexes_o1 <- density_indexes_o1[
    (num_nodes_o1 >= min_ego_nodes) & (num_edges_o1 >= min_ego_edges)
  ]
  density_indexes_o2 <- density_indexes_o2[
    (num_nodes_o2 >= min_ego_nodes) & (num_edges_o2 >= min_ego_edges)
    ]
  # Average manually verified scaled reference counts across density bins
  density_binned_reference_counts_o1 <- rbind(
    (scaled_reference_counts_o1[1,] + scaled_reference_counts_o1[2,]) / 2,
    (scaled_reference_counts_o1[4,] + scaled_reference_counts_o1[5,] + 
      scaled_reference_counts_o1[6,]) / 3,
    ( scaled_reference_counts_o1[3,] +
       scaled_reference_counts_o1[7,] + 
       scaled_reference_counts_o1[8,]) / 3
  )
  rownames(density_binned_reference_counts_o1) <- 1:3
  density_binned_reference_counts_o2 <- rbind(
    (scaled_reference_counts_o2[1,] + scaled_reference_counts_o2[4,]) / 2,
    (scaled_reference_counts_o2[2,] + scaled_reference_counts_o2[6,] +
       scaled_reference_counts_o2[7,]) / 3,
    (scaled_reference_counts_o2[3,] + scaled_reference_counts_o2[5,]) / 2,
    (scaled_reference_counts_o2[8,] + scaled_reference_counts_o2[9,] + 
       scaled_reference_counts_o2[10,]) / 3
  )
  rownames(density_binned_reference_counts_o2) <- 1:4
  
  # Helper functions to calculate expected expected graphlet counts
  expected_expected_graphlet_counts_o1_fn <- function(density_index, node_count) {
    reference_counts <- density_binned_reference_counts_o1[density_index,]
    reference_counts * choose(node_count, graphlet_sizes)
  }
  expected_expected_graphlet_counts_o2_fn <- function(density_index, node_count) {
    reference_counts <- density_binned_reference_counts_o2[density_index,]
    reference_counts * choose(node_count, graphlet_sizes)
  }
  # Calculate expected graphlet counts
  expected_expected_graphlet_counts_ego_o1 <- t(simplify2array(purrr::map2(
    density_indexes_o1, num_nodes_o1[(num_nodes_o1 >= min_ego_nodes)],
    expected_expected_graphlet_counts_o1_fn
  )))
  rownames(expected_expected_graphlet_counts_ego_o1) <- 
    names(ego_networks_o1[(num_nodes_o1 >= min_ego_nodes)])
  expected_expected_graphlet_counts_ego_o2 <- t(simplify2array(purrr::map2(
    density_indexes_o2, num_nodes_o2[(num_nodes_o2 >= min_ego_nodes)],
    expected_expected_graphlet_counts_o2_fn
  )))
  rownames(expected_expected_graphlet_counts_ego_o2) <- 
    names(ego_networks_o2[(num_nodes_o2 >= min_ego_nodes)])
  
  # Sanity check manually derived expected expected counts by comparing against 
  # pre-tested fully applied expected_graphlet_counts_ego function
  expect_equal(expected_expected_graphlet_counts_ego_o1, 
               netdis_expected_graphlet_counts_ego(
                 graph, max_graphlet_size = max_graphlet_size,
                 neighbourhood_size = 1,
                 density_breaks = breaks_o1, 
                 density_binned_reference_counts_o1,
                 min_ego_nodes = min_ego_nodes, min_ego_edges = min_ego_edges
               ))
  expect_equal(expected_expected_graphlet_counts_ego_o2, 
               netdis_expected_graphlet_counts_ego(
                 graph, max_graphlet_size = max_graphlet_size,
                 neighbourhood_size = 2,
                 density_breaks = breaks_o2, 
                 density_binned_reference_counts_o2,
                 min_ego_nodes = min_ego_nodes, min_ego_edges = min_ego_edges
               ))
  
  # Generate partially applied functions using function under test
  actual_expected_graphlet_counts_ego_fn_o1 <- 
    netdis_expected_graphlet_counts_ego_fn(
      graph, max_graphlet_size = max_graphlet_size, neighbourhood_size = 1, 
      min_bin_count = min_bin_count, num_bins = num_bins)
  actual_expected_graphlet_counts_ego_fn_o2 <- 
    netdis_expected_graphlet_counts_ego_fn(
      graph, max_graphlet_size = max_graphlet_size, neighbourhood_size = 2, 
      min_bin_count = min_bin_count, num_bins = num_bins)
  # Generate actual expected accounts by applying generated functions to test
  # graph
  actual_expected_graphlet_counts_ego_o1 <- 
    actual_expected_graphlet_counts_ego_fn_o1(graph)
  actual_expected_graphlet_counts_ego_o2 <- 
    actual_expected_graphlet_counts_ego_fn_o2(graph)
  
  # Compare actual to expected
  expect_equal(actual_expected_graphlet_counts_ego_o1, 
               expected_expected_graphlet_counts_ego_o1)
  expect_equal(actual_expected_graphlet_counts_ego_o2, 
               expected_expected_graphlet_counts_ego_o2)
})
