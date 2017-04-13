
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
  expected_ego_networks_order_1 <- igraph::make_ego_graph(graph, order = 1)
  names(expected_ego_networks_order_1) <- igraph::V(graph)$name
  expected_ego_networks_order_2 <- igraph::make_ego_graph(graph, order = 2)
  names(expected_ego_networks_order_2) <- igraph::V(graph)$name
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
  
  # TODO: Decide if we will generate gdd from scaled ego-network counts. If we do,
  # alter the test below to test this
  # # Test that gdd method gives the expected graphlet degree distributions
  # # 1-step ego-networks
  # actual_gdd_order_1 <- gdd(graph, feature_type = "graphlet", 
  #                           max_graphlet_size = 4, ego_neighbourhood_size = 1)
  # expected_gdd_order_1 <- list(
  #   G0 = dhist(locations = c(1, 4, 5, 6, 7), masses = c(2, 1, 2, 3, 2)),
  #   G1 = dhist(locations = c(0, 2, 3, 5), masses = c(4, 2, 2, 2)),
  #   G2 = dhist(locations = c(0, 1, 2, 4), masses = c(2, 2, 2, 4)),
  #   G3 = dhist(locations = c(0), masses = c(10)),
  #   G4 = dhist(locations = c(0, 1, 2), masses = c(8, 1, 1)),
  #   G5 = dhist(locations = c(0), masses = c(10)),
  #   G6 = dhist(locations = c(0, 1, 2, 3), masses = c(5, 1, 2, 2)),
  #   G7 = dhist(locations = c(0, 1), masses = c(8, 2)),
  #   G8 = dhist(locations = c(0, 1), masses = c(6, 4))
  # )
  # expect_equal(actual_gdd_order_1, expected_gdd_order_1)
  # # 2-step ego-networks
  # actual_gdd_order_2 <- gdd(graph, feature_type = "graphlet", 
  #                           max_graphlet_size = 4, ego_neighbourhood_size = 2)
  # expected_gdd_order_2 <- list(
  #   G0 = dhist(locations = c(5, 8, 9, 10, 11, 13, 15), masses = c(2, 1, 2, 1, 1, 2, 1)),
  #   G1 = dhist(locations = c(5, 8, 10, 13, 14, 18), masses = c(2, 2, 2, 2, 1, 1)),
  #   G2 = dhist(locations = c(1, 2, 4, 5, 6), masses = c(2, 2, 2, 1, 3)),
  #   G3 = dhist(locations = c(0, 4, 6, 10, 11, 15, 21), masses = c(2, 2, 1, 1, 1, 2, 1)),
  #   G4 = dhist(locations = c(0, 1, 2, 3), masses = c(3, 2, 2, 3)),
  #   G5 = dhist(locations = c(0, 1), masses = c(3, 7)),
  #   G6 = dhist(locations = c(2, 4, 5, 6, 8, 9, 11), masses = c(2, 1, 1, 2, 1, 2, 1)),
  #   G7 = dhist(locations = c(0, 1), masses = c(5, 5)),
  #   G8 = dhist(locations = c(0, 1), masses = c(4, 6))
  # )
  # expect_equal(actual_gdd_order_2, expected_gdd_order_2)
})
