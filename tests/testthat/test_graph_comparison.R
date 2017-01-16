library("purrr")
context("Graph comparison")

# GRAPH_CROSS_COMPARISON_SPEC
test_that("Correct cross-comparison specification is generated for virus PPI data", {
  # Load viurs PPI network data in ORCA-compatible edge list format
  data("virusppi")
  
  expected_name_A <- c(rep("EBV", 4), rep("ECL", 3), rep("HSV", 2), 
                       rep("KSHV", 1), rep("VZV", 0))
  expected_index_A <- c(rep(1, 4), rep(2, 3), rep(3, 2), rep(4, 1), rep(5, 0))
  expected_name_B <- c(c("ECL", "HSV", "KSHV", "VZV"), c("HSV", "KSHV", "VZV"), 
                       c("KSHV", "VZV"), c("VZV"))
  expected_index_B <- c(c(2, 3, 4, 5), c(3, 4, 5), c(4, 5), c(5))
  expected <- as.data.frame(cbind(expected_name_A, expected_name_B, 
                             expected_index_A, expected_index_B))
  colnames(expected) <- c("name_a", "name_b", "index_a", "index_b")

  actual <- graph_cross_comparison_spec(virusppi)
  
  matched_output <- function(actual, expected) {
    dims_match <- all.equal(dim(as.matrix(expected)), dim(as.matrix(actual)))
    data_matches <- all.equal(as.matrix(expected), as.matrix(actual))
    headers_match <- all.equal(colnames(expected), colnames(actual))
    return(dims_match && data_matches && headers_match)
  }

  # Check that actual output matches one of the two acceptable outputs at each 
  # cell
  expect_true(matched_output(actual, expected))
})


# # EMD_LP and EMD_CS: Real data tests
# test_that("emd_cs and emd_lp give same output when comapring virus PPI graphs", {
#   # Load viurs PPI network data in ORCA-compatible edge list format
#   data("virusppi")
#   data_indexes <- 1:length(virusppi)
#   data_names <- attr(virusppi, "name")
#   
#   # Calculate ORCA graphlet orbit degree distributions up to graphlet order 4
#   orb_counts <- purrr::map(virusppi, orca::count4)
#   
#   # Set up cross comparison
#   
#   
#   net_emd_for_pair <- function(index_pair, graph_edges) {
#     net_emd()
#   }
#   graph_combinations <- expand.grid(data_indexes)
#   
#   net_emds_lp <- purrr::map(data_indexes)
#   expect_true(FALSE)
# })