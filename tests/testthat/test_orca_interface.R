data_dir <- system.file(file.path("extdata", "VRPINS"), package = "netdist")

context("ORCA interface: Graph to indexed edge round trip")
test_that("Graph to indexed edge list round trip conversion works", {
  g_orig <- igraph::read_graph(file = file.path(data_dir, "EBV-1.txt"), format = "ncol")
  g_rtrip <- netdist::indexed_edges_to_graph(graph_to_indexed_edges(g_orig))
  expect_true(all.equal(igraph::get.edgelist(g_orig),igraph::get.edgelist(g_orig)))
})

context("ORCA interface: Orbit to graphlet counts")
test_that("orbit_to_graphlet_counts summation works", {
  graph <- netdist::virusppi$EBV
  edges <- graph_to_indexed_edges(graph)
  orbit_counts_4 <- orca::count4(edges)
  orbit_counts_5 <- orca::count5(edges)
  # Define orbit indexes belonging to each graphlet using the xero-based indexing
  # from the journal papers, adding one to conver tot he one-based indexing of R
  g0_indexes <- c(0) + 1
  g1_indexes <- c(1:2) + 1
  g2_indexes <- c(3) + 1
  g3_indexes <- c(4:5) + 1
  g4_indexes <- c(6:7) + 1
  g5_indexes <- c(8) + 1
  g6_indexes <- c(9:11) + 1
  g7_indexes <- c(12:13) + 1
  g8_indexes <- c(14) + 1
  g9_indexes <- c(15:17) + 1
  g10_indexes <- c(18:21) + 1
  g11_indexes <- c(22:23) + 1
  g12_indexes <- c(24:26) + 1
  g13_indexes <- c(27:30) + 1
  g14_indexes <- c(31:33) + 1
  g15_indexes <- c(34) + 1
  g16_indexes <- c(35:38) + 1
  g17_indexes <- c(39:42) + 1
  g18_indexes <- c(43:44) + 1
  g19_indexes <- c(45:48) +1
  g20_indexes <- c(49:50) + 1
  g21_indexes <- c(51:53) + 1
  g22_indexes <- c(54:55) + 1
  g23_indexes <- c(56:58) + 1
  g24_indexes <- c(59:61) + 1
  g25_indexes <- c(62:64) + 1
  g26_indexes <- c(65:67) + 1
  g27_indexes <- c(68:69) + 1
  g28_indexes <- c(70:71) + 1
  g29_indexes <- c(72) + 1
  # Get counts for each graphlet
  g0_counts <- rowSums(orbit_counts_5[, g0_indexes, drop = FALSE])
  g1_counts <- rowSums(orbit_counts_5[, g1_indexes, drop = FALSE])
  g2_counts <- rowSums(orbit_counts_5[, g2_indexes, drop = FALSE])
  g3_counts <- rowSums(orbit_counts_5[, g3_indexes, drop = FALSE])
  g4_counts <- rowSums(orbit_counts_5[, g4_indexes, drop = FALSE])
  g5_counts <- rowSums(orbit_counts_5[, g5_indexes, drop = FALSE])
  g6_counts <- rowSums(orbit_counts_5[, g6_indexes, drop = FALSE])
  g7_counts <- rowSums(orbit_counts_5[, g7_indexes, drop = FALSE])
  g8_counts <- rowSums(orbit_counts_5[, g8_indexes, drop = FALSE])
  g9_counts <- rowSums(orbit_counts_5[, g9_indexes, drop = FALSE])
  g10_counts <- rowSums(orbit_counts_5[, g10_indexes, drop = FALSE])
  g11_counts <- rowSums(orbit_counts_5[, g11_indexes, drop = FALSE])
  g12_counts <- rowSums(orbit_counts_5[, g12_indexes, drop = FALSE])
  g13_counts <- rowSums(orbit_counts_5[, g13_indexes, drop = FALSE])
  g14_counts <- rowSums(orbit_counts_5[, g14_indexes, drop = FALSE])
  g15_counts <- rowSums(orbit_counts_5[, g15_indexes, drop = FALSE])
  g16_counts <- rowSums(orbit_counts_5[, g16_indexes, drop = FALSE])
  g17_counts <- rowSums(orbit_counts_5[, g17_indexes, drop = FALSE])
  g18_counts <- rowSums(orbit_counts_5[, g18_indexes, drop = FALSE])
  g19_counts <- rowSums(orbit_counts_5[, g19_indexes, drop = FALSE])
  g20_counts <- rowSums(orbit_counts_5[, g20_indexes, drop = FALSE])
  g21_counts <- rowSums(orbit_counts_5[, g21_indexes, drop = FALSE])
  g22_counts <- rowSums(orbit_counts_5[, g22_indexes, drop = FALSE])
  g23_counts <- rowSums(orbit_counts_5[, g23_indexes, drop = FALSE])
  g24_counts <- rowSums(orbit_counts_5[, g24_indexes, drop = FALSE])
  g25_counts <- rowSums(orbit_counts_5[, g25_indexes, drop = FALSE])
  g26_counts <- rowSums(orbit_counts_5[, g26_indexes, drop = FALSE])
  g27_counts <- rowSums(orbit_counts_5[, g27_indexes, drop = FALSE])
  g28_counts <- rowSums(orbit_counts_5[, g28_indexes, drop = FALSE])
  g29_counts <- rowSums(orbit_counts_5[, g29_indexes, drop = FALSE])
  # Define expected graphlet count matrix for graphlets up to 5 nodes
  expected_graphlet_counts_5 <- 
    cbind(g0_counts, g1_counts, g2_counts, g3_counts, g4_counts, g5_counts, 
          g6_counts, g7_counts, g8_counts, g9_counts, g10_counts, g11_counts,
          g12_counts, g13_counts, g14_counts, g15_counts, g16_counts,
          g17_counts, g18_counts, g19_counts, g20_counts, g21_counts,
          g22_counts, g23_counts, g24_counts, g25_counts, g26_counts,
          g27_counts, g28_counts, g29_counts)
  colnames(expected_graphlet_counts_5) <- 
             c("G0", "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9","G10",
               "G11", "G12", "G13", "G14", "G15", "G16", "G17", "G18", "G19", 
               "G20", "G21", "G22", "G23", "G24", "G25", "G26", "G27", "G28", 
               "G29")
  # Define epected graphlet count matrix for graphlets up to 4 nodes by selecting
  # a subset of the matrix for graphlets up to 5 nodes
  expected_graphlet_counts_4 <- expected_graphlet_counts_5[,1:9]
  # Calculate actual graphlet counts from functions under test
  actual_graphlet_counts_4 <- orbit_to_graphlet_counts(orbit_counts_4)
  actual_graphlet_counts_5 <- orbit_to_graphlet_counts(orbit_counts_5)
  # Check expected and actual graphlet counts match
  expect_equal(actual_graphlet_counts_4, expected_graphlet_counts_4)
  expect_equal(actual_graphlet_counts_5, expected_graphlet_counts_5)
})


context("ORCA interface: Graphlet-based degree distributions")
test_that("gdd works", {
  graph <- netdist::virusppi$EBV
  edges <- graph_to_indexed_edges(graph)
  # Caclulate expected outputs (NOTE: relies on orbit_to_graphlet_counts and
  # orca_counts_to_graphlet_orbit_degree_distribution methods)
  orbit_counts_4 <- orca::count4(edges)
  orbit_counts_5 <- orca::count5(edges)
  graphlet_counts_4 <- orbit_to_graphlet_counts(orbit_counts_4)
  graphlet_counts_5 <- orbit_to_graphlet_counts(orbit_counts_5)
  gdd_orbit_4_expected <- orca_counts_to_graphlet_orbit_degree_distribution(orbit_counts_4)
  gdd_orbit_5_expected <- orca_counts_to_graphlet_orbit_degree_distribution(orbit_counts_5)
  gdd_graphlet_4_expected <- orca_counts_to_graphlet_orbit_degree_distribution(graphlet_counts_4)
  gdd_graphlet_5_expected <- orca_counts_to_graphlet_orbit_degree_distribution(graphlet_counts_5)
  # Calculate actual outputs from code under test
  gdd_orbit_4_actual <- gdd(graph, feature_type = "orbit", max_graphlet_size = 4)
  gdd_orbit_5_actual <- gdd(graph, feature_type = "orbit", max_graphlet_size = 5)
  gdd_graphlet_4_actual <- gdd(graph, feature_type = "graphlet", max_graphlet_size = 4)
  gdd_graphlet_5_actual <- gdd(graph, feature_type = "graphlet", max_graphlet_size = 5)
  gdd_default_4_actual <- gdd(graph, max_graphlet_size = 4)
  gdd_default_5_actual <- gdd(graph, max_graphlet_size = 5)
  gdd_orbit_default_actual <- gdd(graph, feature_type = "orbit")
  gdd_graphlet_default_actual <- gdd(graph, feature_type = "graphlet")
  gdd_default_default_actual <- gdd(graph)
  # Compare actual with expected
  expect_equal(gdd_orbit_4_actual, gdd_orbit_4_expected)
  expect_equal(gdd_orbit_5_actual, gdd_orbit_5_expected)
  expect_equal(gdd_graphlet_4_actual, gdd_graphlet_4_expected)
  expect_equal(gdd_graphlet_5_actual, gdd_graphlet_5_expected)
  expect_equal(gdd_default_4_actual, gdd_orbit_4_expected)
  expect_equal(gdd_default_5_actual, gdd_orbit_5_expected)
  expect_equal(gdd_orbit_default_actual, gdd_orbit_4_expected)
  expect_equal(gdd_graphlet_default_actual, gdd_graphlet_4_expected)
  expect_equal(gdd_default_default_actual, gdd_orbit_4_expected)
})

context("ORCA interface: Ego-network graphlet outputs for manually verified networks")
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
    count_graphlets_ego(graph, max_graphlet_size = max_graphlet_size, 
                        neighbourhood_size = 1)
  actual_counts_order_2 <- 
    count_graphlets_ego(graph, max_graphlet_size = max_graphlet_size, 
                        neighbourhood_size = 2)
  
  # Set manually verified counts
  # 1-step ego networks
  expected_counts_order_1 <- rbind(
    c(6, 5, 2, 0, 1, 0, 2, 1, 0),
    c(5, 5, 1, 0, 2, 0, 2, 0, 0),
    c(1, 0, 0, 0, 0, 0, 0, 0, 0),
    c(5, 2, 2, 0, 0, 0, 0, 1, 0),
    c(1, 0, 0, 0, 0, 0, 0, 0, 0),
    c(4, 2, 1, 0, 0, 0, 1, 0, 0),
    c(7, 3, 4, 0, 0, 0, 3, 0, 1),
    c(7, 3, 4, 0, 0, 0, 3, 0, 1),
    c(6, 0, 4, 0, 0, 0, 0, 0, 1),
    c(6, 0, 4, 0, 0, 0, 0, 0, 1)
  )
  rownames(expected_counts_order_1) <- node_labels
  colnames(expected_counts_order_1) <- graphlet_labels
  # 2-step ego networks
  expected_counts_order_2 <- rbind(
    c(15, 18, 6, 21, 3, 1, 11, 1, 1),
    c( 8, 10, 2,  6, 3, 0,  4, 1, 0),
    c( 5,  5, 1,  0, 2, 0,  2, 0, 0),
    c(10, 14, 2, 11, 3, 1,  5, 1, 0),
    c( 5,  5, 1,  0, 2, 0,  2, 0, 0),
    c(13, 13, 6, 15, 1, 1,  9, 1, 1),
    c(13, 13, 6, 15, 1, 1,  9, 1, 1),
    c(11, 10, 5, 10 ,0 ,1,  8, 0, 1),
    c( 9,  8, 4,  4, 0, 1,  6, 0, 1),
    c( 9,  8, 4,  4, 0, 1,  6, 0, 1)
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
    count_graphlets_ego(graph, max_graphlet_size = max_graphlet_size, 
                        neighbourhood_size = 1, return_ego_networks = TRUE)
  actual_counts_with_networks_order_2 <- 
    count_graphlets_ego(graph, max_graphlet_size = max_graphlet_size, 
                        neighbourhood_size = 2, return_ego_networks = TRUE)
  # 3. Compare
  expect_equal(actual_counts_with_networks_order_1, expected_counts_with_networks_order_1)
  expect_equal(actual_counts_with_networks_order_2, expected_counts_with_networks_order_2)
  
  # Test that gdd method gives the expected graphlet degree distributions
  # 1-step ego-networks
  actual_gdd_order_1 <- gdd(graph, feature_type = "graphlet", 
                            max_graphlet_size = 4, ego_neighbourhood_size = 1)
  expected_gdd_order_1 <- list(
    G0 = dhist(locations = c(1, 4, 5, 6, 7), masses = c(2, 1, 2, 3, 2)),
    G1 = dhist(locations = c(0, 2, 3, 5), masses = c(4, 2, 2, 2)),
    G2 = dhist(locations = c(0, 1, 2, 4), masses = c(2, 2, 2, 4)),
    G3 = dhist(locations = c(0), masses = c(10)),
    G4 = dhist(locations = c(0, 1, 2), masses = c(8, 1, 1)),
    G5 = dhist(locations = c(0), masses = c(10)),
    G6 = dhist(locations = c(0, 1, 2, 3), masses = c(5, 1, 2, 2)),
    G7 = dhist(locations = c(0, 1), masses = c(8, 2)),
    G8 = dhist(locations = c(0, 1), masses = c(6, 4))
  )
  expect_equal(actual_gdd_order_1, expected_gdd_order_1)
  # 2-step ego-networks
  actual_gdd_order_2 <- gdd(graph, feature_type = "graphlet", 
                            max_graphlet_size = 4, ego_neighbourhood_size = 2)
  expected_gdd_order_2 <- list(
    G0 = dhist(locations = c(5, 8, 9, 10, 11, 13, 15), masses = c(2, 1, 2, 1, 1, 2, 1)),
    G1 = dhist(locations = c(5, 8, 10, 13, 14, 18), masses = c(2, 2, 2, 2, 1, 1)),
    G2 = dhist(locations = c(1, 2, 4, 5, 6), masses = c(2, 2, 2, 1, 3)),
    G3 = dhist(locations = c(0, 4, 6, 10, 11, 15, 21), masses = c(2, 2, 1, 1, 1, 2, 1)),
    G4 = dhist(locations = c(0, 1, 2, 3), masses = c(3, 2, 2, 3)),
    G5 = dhist(locations = c(0, 1), masses = c(3, 7)),
    G6 = dhist(locations = c(2, 4, 5, 6, 8, 9, 11), masses = c(2, 1, 1, 2, 1, 2, 1)),
    G7 = dhist(locations = c(0, 1), masses = c(5, 5)),
    G8 = dhist(locations = c(0, 1), masses = c(4, 6))
  )
  expect_equal(actual_gdd_order_2, expected_gdd_order_2)
})
