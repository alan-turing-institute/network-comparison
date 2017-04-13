
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
  graph_n11_ego1 <- igraph::make_ego_graph(graph_n11, order = 1)
  graph_n37_ego1 <- igraph::make_ego_graph(graph_n37, order = 1)
  graph_n73_ego1 <- igraph::make_ego_graph(graph_n73, order = 1)
  graph_n11_ego2 <- igraph::make_ego_graph(graph_n11, order = 2)
  graph_n37_ego2 <- igraph::make_ego_graph(graph_n37, order = 2)
  graph_n73_ego2 <- igraph::make_ego_graph(graph_n73, order = 2)
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
