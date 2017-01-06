library("netdist")
library("purrr")
context("NetEMD")


# BIN NORMALISATION: Property based tests
test_that("normalise_histogram_mass output sums to 1", {
  # Generate histograms with random masses and random centres (uniformly random)
  num_hists <- 10
  num_bins <- 100
  
  mass_min <- 0
  mass_max <- 100
  rand_bin_masses <- function() {return(runif(num_bins, mass_min, mass_max))}
  bin_mass_lists <- replicate(num_hists, rand_bin_masses(), simplify = FALSE)
  
  actual_sums <- purrr::map_dbl(bin_mass_lists, function(bin_masses) {sum(normalise_histogram_mass(bin_masses))})
  expected <- 1
  purrr::map_dbl(actual_sums, function(actual) {expect_equal(actual, expected)})
})

test_that("normalise_histogram_variance output has variance of 1", {
  expect_true(FALSE)
})

# COST_MATRIX: Property-based tests
test_that("cost_matrix returns all zeros when all bin locations are identical", {
  bin_centres1 <- c(1, 1, 1, 1, 1, 1, 1)
  bin_centres2 <- bin_centres1
  expected <- matrix(0, nrow = length(bin_centres1), ncol = length(bin_centres2))
  expect_equal(cost_matrix(bin_centres1, bin_centres2), expected)
})

test_that("cost_matrix returns zeros along diagonal when both sets of bin locations are the same", {
  bin_centres1 <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)
  bin_centres2 <- bin_centres1
  expected <- rep(0, length(bin_centres1))
  expect_equal(diag(cost_matrix(bin_centres1, bin_centres2)), expected)
})

test_that("cost_matrix returns zeros along diagonal and taxicab distance from all zeros for all other elements when both sets of bin locations are the same and are a sequence of consecutive integers", {
  bin_centres1 <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)
  bin_centres2 <- bin_centres1
  num_bins <- length(bin_centres1)
  expected <- toeplitz(1:num_bins)-1
  expect_equal(cost_matrix(bin_centres1, bin_centres2), expected)
})

test_that("cost_matrix is correct size when the two histograms are of different lengths", {
  bin_centres1 <- c(1, 2, 3, 4, 5, 6, 7)
  bin_centres2 <- c(8, 9, 10)
  
  cm <- cost_matrix(bin_centres1, bin_centres2)
  
  expect_equal(nrow(cm), length(bin_centres1))
  expect_equal(ncol(cm), length(bin_centres2))
})

# AUGMENT_HISTOGRAMS
test_that("augment_histograms works A", {
  bin_masses1 <- c(1, 1, 1)
  bin_masses2 <- c(1, 1, 1)
  bin_centres1 <- c(1, 3, 5)
  bin_centres2 <- c(2, 4, 6)
  
  expected <- list(
    bin_masses1 = c(1, 1, 1, 0, 0, 0),
    bin_masses2 = c(1, 1, 1, 0, 0, 0),
    bin_centres1 = c(1, 3, 5, 2, 4, 6),
    bin_centres2 = c(2, 4, 6, 1, 3, 5)
  )
  actual <- augment_histograms(bin_masses1, bin_masses2, bin_centres1, bin_centres2)
  expect_equal(actual, expected)
  
})
test_that("augment_histograms works B", {
  bin_masses1 <- c(1, 1, 1)
  bin_masses2 <- c(1, 1, 1)
  bin_centres1 <- c(1, 3, 5)
  bin_centres2 <- c(4, 5, 6)
  
  expected <- list(
    bin_masses1 = c(1, 1, 1, 0, 0),
    bin_masses2 = c(1, 1, 1, 0, 0),
    bin_centres1 = c(1, 3, 5, 4, 6),
    bin_centres2 = c(4, 5, 6, 1, 3)
  )
  actual <- augment_histograms(bin_masses1, bin_masses2, bin_centres1, bin_centres2)
  expect_equal(actual, expected)
})

# EMD_LP and EMD_CS: Property-based tests
test_that("EMD methods return 0 when comparing a 1D feature distribution to itself",{
  bin_masses1 <- c(0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0)
  bin_masses2 <- bin_masses1
  bin_centres1 <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)
  bin_centres2 <- bin_centres1
  
  expected <- 0
  expect_equal(emd_lp(bin_masses1, bin_masses2, bin_centres1, bin_centres2), expected)
  expect_equal(emd_cs(bin_masses1, bin_masses2, bin_centres1, bin_centres2), expected)
})

test_that("EMD methods return numBins/2 when offsetting a symmetric discrete triangle distribution by 1", {
  cost_fn <- function(triangle_width) {
    move_dist <- ceiling((triangle_width+1)/2)
    num_moves <- ceiling(triangle_width/2)
    return(move_dist * num_moves)
  }
  
  bin_masses1 <- c(0, 1, 2, 3, 4, 4, 3, 2, 1, 0)
  bin_masses2 <- c(0, 0, 1, 2, 3, 4, 4, 3, 2, 1)
  bin_centres1 <- 1:length(bin_masses1)
  bin_centres2 <- 1:length(bin_masses2)
  num_nonzero_bins <- sum(bin_masses1 > 0)
  expected <- cost_fn(num_nonzero_bins)
  expect_equal(emd_lp(bin_masses1, bin_masses2, bin_centres1, bin_centres2), expected)
  expect_equal(emd_cs(bin_masses1, bin_masses2, bin_centres1, bin_centres2), expected)
  
  bin_masses1 <- c(0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0)
  bin_masses2 <- c(0, 0, 1, 2, 3, 4, 5, 4, 3, 2, 1)
  bin_centres1 <- 1:length(bin_masses1)
  bin_centres2 <- 1:length(bin_masses2)
  expected <- cost_fn(sum(bin_masses1 > 0))
  expect_equal(emd_lp(bin_masses1, bin_masses2, bin_centres1, bin_centres2), expected)
  expect_equal(emd_cs(bin_masses1, bin_masses2, bin_centres1, bin_centres2), expected)
  
  bin_masses1 <- c(0, 1, 2, 3, 4, 5, 5, 4, 3, 2, 1, 0)
  bin_masses2 <- c(0, 0, 1, 2, 3, 4, 5, 5, 4, 3, 2, 1)
  bin_centres1 <- 1:length(bin_masses1)
  bin_centres2 <- 1:length(bin_masses2)
  expected <- cost_fn(sum(bin_masses1 > 0))
  expect_equal(emd_lp(bin_masses1, bin_masses2, bin_centres1, bin_centres2), expected)
  expect_equal(emd_cs(bin_masses1, bin_masses2, bin_centres1, bin_centres2), expected)
  
  bin_masses1 <- c(0, 1, 2, 3, 4, 5, 6, 5, 4, 3, 2, 1, 0)
  bin_masses2 <- c(0, 0, 1, 2, 3, 4, 5, 6, 5, 4, 3, 2, 1)
  bin_centres1 <- 1:length(bin_masses1)
  bin_centres2 <- 1:length(bin_masses2)
  expected <- cost_fn(sum(bin_masses1 > 0))
  expect_equal(emd_lp(bin_masses1, bin_masses2, bin_centres1, bin_centres2), expected)
  expect_equal(emd_cs(bin_masses1, bin_masses2, bin_centres1, bin_centres2), expected)
  
})

test_that("EMD methods return same result for densely and sparsely specified bins", {
  sparse_bin_masses1 <- c(1, 1, 1, 1, 1, 1)
  sparse_bin_masses2 <- c(1, 1, 1, 1, 1, 1)
  sparse_bin_centres1 <- c(1, 2, 4, 7, 11, 16)
  sparse_bin_centres2 <- c(21, 22, 24, 27, 31, 36)
  
  dense_bin_centres1 <- 1:36
  dense_bin_centres2 <- dense_bin_centres1
  bin_mass_sequence <- c(1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1)
  bin_mass_padding <- rep(0, length(dense_bin_centres1) - length(bin_mass_sequence))
  dense_bin_masses1 <- c(bin_mass_sequence, bin_mass_padding)
  dense_bin_masses2 <- c(bin_mass_padding, bin_mass_sequence)
  
  expect_equal(emd_lp(dense_bin_masses1, dense_bin_masses2, 
                    dense_bin_centres1, dense_bin_centres2),
              emd_lp(sparse_bin_masses1, sparse_bin_masses2, 
                    sparse_bin_centres1, sparse_bin_centres2))
  expect_equal(emd_cs(dense_bin_masses1, dense_bin_masses2, 
                    dense_bin_centres1, dense_bin_centres2),
              emd_cs(sparse_bin_masses1, sparse_bin_masses2, 
                    sparse_bin_centres1, sparse_bin_centres2))
})

test_that("EMD methods return same result when order of densely specified bins is changed", {
  bin_masses1 <- c(1, 1, 1, 1, 0, 0, 0, 0, 0)
  bin_masses2 <- c(0, 0, 0, 0, 0, 1, 1, 1, 1)
  bin_centres1 <- 1:length(bin_masses1)
  bin_centres2 <- 1:length(bin_masses2)
  
  permuted_indexes1 <- sample(1:length(bin_masses1), replace = FALSE)
  permuted_indexes2 <- sample(1:length(bin_masses2), replace = FALSE)
  
  permuted_bin_masses1 <- bin_masses1[permuted_indexes1]
  permuted_bin_centres1 <- bin_centres1[permuted_indexes1]
  permuted_bin_masses2 <- bin_masses2[permuted_indexes2]
  permuted_bin_centres2 <- bin_centres2[permuted_indexes2]
  
  expect_equal(emd_lp(bin_masses1, bin_masses2, bin_centres1, bin_centres2),
              emd_lp(permuted_bin_masses1, permuted_bin_masses2, 
                    permuted_bin_centres1, permuted_bin_centres2))
  expect_equal(emd_cs(bin_masses1, bin_masses2, bin_centres1, bin_centres2),
              emd_cs(permuted_bin_masses1, permuted_bin_masses2, 
                    permuted_bin_centres1, permuted_bin_centres2))
})

test_that("EMD methods return same result when order of sparsely specified bins is changed", {
  bin_masses1 <- c(1, 1, 1, 1, 1, 1)
  bin_masses2 <- c(1, 1, 1, 1, 1, 1)
  bin_centres1 <- c(1, 2, 4, 8, 16, 32)
  bin_centres2 <- c(-32, -16, -8, -4, -2, -1)
  
  permuted_indexes1 <- sample(1:length(bin_masses1), replace = FALSE)
  permuted_indexes2 <- sample(1:length(bin_masses2), replace = FALSE)
  
  permuted_bin_masses1 <- bin_masses1[permuted_indexes1]
  permuted_bin_centres1 <- bin_centres1[permuted_indexes1]
  permuted_bin_masses2 <- bin_masses2[permuted_indexes2]
  permuted_bin_centres2 <- bin_centres2[permuted_indexes2]
  
  expect_equal(emd_lp(bin_masses1, bin_masses2, bin_centres1, bin_centres2),
               emd_lp(permuted_bin_masses1, permuted_bin_masses2, 
                      permuted_bin_centres1, permuted_bin_centres2))
  expect_equal(emd_cs(bin_masses1, bin_masses2, bin_centres1, bin_centres2),
               emd_cs(permuted_bin_masses1, permuted_bin_masses2, 
                      permuted_bin_centres1, permuted_bin_centres2))
  })

# NetEMD: Property-based tests
test_that("net_emd returns 0 when comparing any histogram offset against itself", {
  expect_true(FALSE)
})

# EMD_LP and EMD_CS: Real data tests
test_that("EMD methods return correct results for sample virus PPI data sets", {
  file_names <- c("EBV-1.txt", "ECL-1.txt", "HSV-1-1.txt", "KSHV-1.txt", "VZV-1.txt")
  dataset_names <- c("EBV", "ECL", "HSV", "KSHV", "VZV")
  file_formats <- c("ncol", "ncol", "ncol", "ncol", "ncol")
  data <- data.frame(dataset_names, file_names, file_formats)
})

