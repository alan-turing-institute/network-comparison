library("purrr")

# GRAPH_CROSS_COMPARISON_SPEC
test_that("Correct cross-comparison specification is generated for virus PPI data", {
  # Load viurs PPI network data in ORCA-compatible edge list format
  data("virusppi")
  
  expected_name_1 <- c(rep("EBV", 5), rep("ECL", 5), rep("HSV", 5), 
                       rep("KSHV", 5), rep("VZV", 5))
  expected_index_1 <- c(rep(1, 5), rep(2, 5), rep(3, 5), rep(4, 5), rep(5, 5))
  expected_name_2 <- c(rep(c("EBV", "ECL", "HSV", "KSHV", "VZV"), 5))
  expected_index_2 <- c(rep(c(1, 2, 3, 4 ,5), 5))
  
  expected_1_first <- as.data.frame(cbind(expected_name_1, expected_name_2, 
                             expected_index_1, expected_index_2))
  expected_2_first <- as.data.frame(cbind(expected_name_2, expected_name_1, 
                             expected_index_2, expected_index_1))
  colnames(expected_1_first) <- c("Name A", "Name B", "Index A", "Index B")
  colnames(expected_2_first) <- colnames(expected_1_first)
  
  
  actual <- graph_cross_comparison_spec(virusppi)
  
  matched_output <- function(actual, expected) {
    data_matches <- all.equal(as.matrix(expected_2_first), as.matrix(actual))
    headers_match <- all.equal(colnames(expected), colnames(actual))
  }
  
  # Check that actual output matches one of the two acceptable outputs at each 
  # cell
  expect_true(matched_output(actual, expected_1_first) 
              | matched_output(actual, expected_2_first))
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