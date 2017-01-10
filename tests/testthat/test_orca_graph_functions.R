library("igraph")
context("ORCA Graph Functions")

data_dir <- system.file(file.path("inst", "extdata", "VRPINS"), package = "netdist")

test_that("Graph to indexed edge list round trip conversion works", {
  g_orig <- igraph::read_graph(file = file.path(data_dir, "EBV-1.txt"), format = "ncol")
  g_rtrip <- indexed_edges_to_graph(graph_to_indexed_edges(g_orig))
  expect_true(all.equal(igraph::get.edgelist(g_orig),igraph::get.edgelist(g_orig)))
})