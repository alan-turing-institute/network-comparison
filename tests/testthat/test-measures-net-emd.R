library("netdist")
context("NetEMD")

# COST_MATRIX: Property-based tests
test_that("cost_matrix returns all zeros when all bin locations are identical", {
  bin_locations1 <- c(1, 1, 1, 1, 1, 1, 1)
  bin_locations2 <- bin_locations1
  expected <- matrix(0, nrow = length(bin_locations1), ncol = length(bin_locations2))
  expect_equal(cost_matrix(bin_locations1, bin_locations2), expected)
})

test_that("cost_matrix returns zeros along diagonal when both sets of bin locations are the same", {
  bin_locations1 <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)
  bin_locations2 <- bin_locations1
  expected <- rep(0, length(bin_locations1))
  expect_equal(diag(cost_matrix(bin_locations1, bin_locations2)), expected)
})

test_that("cost_matrix returns zeros along diagonal and taxicab distance from all zeros for all other elements when both sets of bin locations are the same and are a sequence of consecutive integers", {
  bin_locations1 <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)
  bin_locations2 <- bin_locations1
  num_bins <- length(bin_locations1)
  expected <- toeplitz(1:num_bins)-1
  expect_equal(cost_matrix(bin_locations1, bin_locations2), expected)
})

# EMD_LP: Property-based tests
test_that("emd_lp returns 0 when comparing a 1D feature distribution to itself (implicit unit interval locations)",{
  bin_counts1 <- c(0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0)
  bin_counts2 <- bin_counts1
  expect_equal(emd_lp(bin_counts1, bin_counts2), 0)
})

test_that("emd_lp returns 0 when comparing a 1D feature distribution to itself (explicit locations)",{
  bin_counts1 <- c(0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0)
  bin_counts2 <- bin_counts1
  bin_locations1 <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)
  bin_locations2 <- bin_locations1
  expect_equal(emd_lp(bin_counts1, bin_counts2, bin_locations1, bin_locations2), 0)
})

test_that("emd_lp returns numBins/2 when offsetting a symmetric discrete triangle distribution by 1 (implicit unit interal locations)", {
  cost_fn <- function(triangle_width) {
    move_dist <- ceiling((triangle_width+1)/2)
    num_moves <- ceiling(triangle_width/2)
    return(move_dist * num_moves)
  }
  
  bin_counts1 <- c(0, 1, 2, 3, 4, 4, 3, 2, 1, 0)
  bin_counts2 <- c(0, 0, 1, 2, 3, 4, 4, 3, 2, 1)
  expected = cost_fn(sum(bin_counts1>0))
  expect_equal(emd_lp(bin_counts1, bin_counts2), expected)
  
  bin_counts1 <- c(0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0)
  bin_counts2 <- c(0, 0, 1, 2, 3, 4, 5, 4, 3, 2, 1)
  expected = cost_fn(sum(bin_counts1>0))
  expect_equal(emd_lp(bin_counts1, bin_counts2), expected)
  
  bin_counts1 <- c(0, 1, 2, 3, 4, 5, 5, 4, 3, 2, 1, 0)
  bin_counts2 <- c(0, 0, 1, 2, 3, 4, 5, 5, 4, 3, 2, 1)
  expected = cost_fn(sum(bin_counts1>0))
  expect_equal(emd_lp(bin_counts1, bin_counts2), expected)
  
  bin_counts1 <- c(0, 1, 2, 3, 4, 5, 6, 5, 4, 3, 2, 1, 0)
  bin_counts2 <- c(0, 0, 1, 2, 3, 4, 5, 6, 5, 4, 3, 2, 1)
  expected = cost_fn(sum(bin_counts1>0))
  expect_equal(emd_lp(bin_counts1, bin_counts2), expected)
  
})
