context("EMD: Cost matrix")
# COST_MATRIX: Property-based tests
test_that("cost_matrix returns all zeros when all bin locations are identical",{
  bin_centres1 <- c(1, 1, 1, 1, 1, 1, 1)
  bin_centres2 <- bin_centres1
  expected <- matrix(0, nrow = length(bin_centres1), 
                     ncol = length(bin_centres2))
  expect_equal(cost_matrix(bin_centres1, bin_centres2), expected)
})

test_that("cost_matrix returns zeros along diagonal when both sets of bin 
          locations are the same", {
            bin_centres1 <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)
            bin_centres2 <- bin_centres1
            expected <- rep(0, length(bin_centres1))
            expect_equal(diag(cost_matrix(bin_centres1, bin_centres2)), expected)
            })

test_that("cost_matrix returns zeros along diagonal and taxicab distance from 
          all zeros for all other elements when both sets of bin locations are  
          the same and are a sequence of consecutive integers", {
            bin_centres1 <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)
            bin_centres2 <- bin_centres1
            num_bins <- length(bin_centres1)
            expected <- toeplitz(1:num_bins)-1
            expect_equal(cost_matrix(bin_centres1, bin_centres2), expected)
            })

test_that("cost_matrix is correct size when the two histograms are of different 
          lengths", {
            bin_centres1 <- c(1, 2, 3, 4, 5, 6, 7)
            bin_centres2 <- c(8, 9, 10)
            
            cm <- cost_matrix(bin_centres1, bin_centres2)
            
            expect_equal(nrow(cm), length(bin_centres1))
            expect_equal(ncol(cm), length(bin_centres2))
            })

context("EMD: EMD")
# EMD: Property-based tests
test_that("EMD methods return 0 when comparing a 1D feature distribution to 
          itself",{
            bin_masses1 <- c(0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0)
            bin_masses2 <- bin_masses1
            bin_centres1 <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)
            bin_centres2 <- bin_centres1
            histogram1 <- dhist(masses = bin_masses1, locations = bin_centres1)
            histogram2 <- dhist(masses = bin_masses2, locations = bin_centres2)
            
            expected <- 0
            expect_equal(emd_lp(bin_masses1, bin_masses2, 
                                bin_centres1, bin_centres2), expected)
            expect_equal(emd_cs(histogram1, histogram2), expected)
            expect_equal(emd(histogram1, histogram2), expected)
            })

test_that("EMD methods return numBins/2 when offsetting a symmetric discrete 
          triangle distribution by 1", {
            cost_fn <- function(triangle_width) {
              move_dist <- ceiling((triangle_width+1)/2)
              num_moves <- ceiling(triangle_width/2)
              return(move_dist * num_moves)
            }
            
            # Triangle(4, even), shifting by changing masses
            bin_masses1 <- c(0, 1, 2, 3, 4, 4, 3, 2, 1, 0)
            bin_masses2 <- c(0, 0, 1, 2, 3, 4, 4, 3, 2, 1)
            bin_centres1 <- 1:length(bin_masses1)
            bin_centres2 <- 1:length(bin_masses2)
            histogram1 <- dhist(masses = bin_masses1, locations = bin_centres1)
            histogram2 <- dhist(masses = bin_masses2, locations = bin_centres2)
            num_nonzero_bins <- sum(bin_masses1 > 0)
            expected <- cost_fn(num_nonzero_bins)
            emd_lp(bin_masses1, bin_masses2, bin_centres1, bin_centres2)
            expect_equal(emd_lp(bin_masses1, bin_masses2, 
                                bin_centres1, bin_centres2), expected)
            expect_equal(emd_cs(histogram1, histogram2), expected)
            expect_equal(emd(histogram1, histogram2), expected)
            
            # Triangle(4, even), shifting by changing centres
            bin_masses1 <- c(0, 1, 2, 3, 4, 4, 3, 2, 1, 0)
            bin_masses2 <- c(0, 1, 2, 3, 4, 4, 3, 2, 1, 0)
            bin_centres1 <- 1:length(bin_masses1)
            bin_centres2 <- 1:length(bin_masses2) + 1
            histogram1 <- dhist(masses = bin_masses1, locations = bin_centres1)
            histogram2 <- dhist(masses = bin_masses2, locations = bin_centres2)
            num_nonzero_bins <- sum(bin_masses1 > 0)
            expected <- cost_fn(num_nonzero_bins)
            expect_equal(emd_lp(bin_masses1, bin_masses2, 
                                bin_centres1, bin_centres2), expected)
            expect_equal(emd_cs(histogram1, histogram2), expected)
            expect_equal(emd(histogram1, histogram2), expected)
            
            # Triangle(5, odd), shifting by changing masses
            bin_masses1 <- c(0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0)
            bin_masses2 <- c(0, 0, 1, 2, 3, 4, 5, 4, 3, 2, 1)
            bin_centres1 <- 1:length(bin_masses1)
            bin_centres2 <- 1:length(bin_masses2)
            histogram1 <- dhist(masses = bin_masses1, locations = bin_centres1)
            histogram2 <- dhist(masses = bin_masses2, locations = bin_centres2)
            expected <- cost_fn(sum(bin_masses1 > 0))
            expect_equal(emd_lp(bin_masses1, bin_masses2, 
                                bin_centres1, bin_centres2), expected)
            expect_equal(emd_cs(histogram1, histogram2), expected)
            expect_equal(emd(histogram1, histogram2), expected)
            
            # Triangle(5, odd), shifting by changing masses
            bin_masses1 <- c(0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0)
            bin_masses2 <- c(0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0)
            bin_centres1 <- 1:length(bin_masses1)
            bin_centres2 <- 1:length(bin_masses2) + 1
            histogram1 <- dhist(masses = bin_masses1, locations = bin_centres1)
            histogram2 <- dhist(masses = bin_masses2, locations = bin_centres2)
            expected <- cost_fn(sum(bin_masses1 > 0))
            expect_equal(emd_lp(bin_masses1, bin_masses2, 
                                bin_centres1, bin_centres2), expected)
            expect_equal(emd_cs(histogram1, histogram2), expected)
            expect_equal(emd(histogram1, histogram2), expected)
            
            # Triangle(5, even), shifting by changing masses
            bin_masses1 <- c(0, 1, 2, 3, 4, 5, 5, 4, 3, 2, 1, 0)
            bin_masses2 <- c(0, 0, 1, 2, 3, 4, 5, 5, 4, 3, 2, 1)
            bin_centres1 <- 1:length(bin_masses1)
            bin_centres2 <- 1:length(bin_masses2)
            histogram1 <- dhist(masses = bin_masses1, locations = bin_centres1)
            histogram2 <- dhist(masses = bin_masses2, locations = bin_centres2)
            expected <- cost_fn(sum(bin_masses1 > 0))
            expect_equal(emd_lp(bin_masses1, bin_masses2, 
                                bin_centres1, bin_centres2), expected)
            expect_equal(emd_cs(histogram1, histogram2), expected)
            expect_equal(emd(histogram1, histogram2), expected)
            
            # Triangle(5, even), shifting by changing centres
            bin_masses1 <- c(0, 1, 2, 3, 4, 5, 5, 4, 3, 2, 1, 0)
            bin_masses2 <- c(0, 1, 2, 3, 4, 5, 5, 4, 3, 2, 1, 0)
            bin_centres1 <- 1:length(bin_masses1)
            bin_centres2 <- 1:length(bin_masses2) + 1
            histogram1 <- dhist(masses = bin_masses1, locations = bin_centres1)
            histogram2 <- dhist(masses = bin_masses2, locations = bin_centres2)
            expected <- cost_fn(sum(bin_masses1 > 0))
            expect_equal(emd_lp(bin_masses1, bin_masses2, 
                                bin_centres1, bin_centres2), expected)
            expect_equal(emd_cs(histogram1, histogram2), expected)
            expect_equal(emd(histogram1, histogram2), expected)
            
            # Triangle(6, odd), shifting by changing masses
            bin_masses1 <- c(0, 1, 2, 3, 4, 5, 6, 5, 4, 3, 2, 1, 0)
            bin_masses2 <- c(0, 0, 1, 2, 3, 4, 5, 6, 5, 4, 3, 2, 1)
            bin_centres1 <- 1:length(bin_masses1)
            bin_centres2 <- 1:length(bin_masses2)
            histogram1 <- dhist(masses = bin_masses1, locations = bin_centres1)
            histogram2 <- dhist(masses = bin_masses2, locations = bin_centres2)
            expected <- cost_fn(sum(bin_masses1 > 0))
            expect_equal(emd_lp(bin_masses1, bin_masses2, 
                                bin_centres1, bin_centres2), expected)
            expect_equal(emd_cs(histogram1, histogram2), expected)
            expect_equal(emd(histogram1, histogram2), expected)
            
            # Triangle(6, odd), shifting by changing centres
            bin_masses1 <- c(0, 1, 2, 3, 4, 5, 6, 5, 4, 3, 2, 1, 0)
            bin_masses2 <- c(0, 1, 2, 3, 4, 5, 6, 5, 4, 3, 2, 1, 0)
            bin_centres1 <- 1:length(bin_masses1)
            bin_centres2 <- 1:length(bin_masses2) + 1
            histogram1 <- dhist(masses = bin_masses1, locations = bin_centres1)
            histogram2 <- dhist(masses = bin_masses2, locations = bin_centres2)
            expected <- cost_fn(sum(bin_masses1 > 0))
            expect_equal(emd_lp(bin_masses1, bin_masses2, 
                                bin_centres1, bin_centres2), expected)
            expect_equal(emd_cs(histogram1, histogram2), expected)
            expect_equal(emd(histogram1, histogram2), expected)
            
            })

test_that("EMD methods return same result for densely and sparsely specified 
          bins", {
            sparse_bin_masses1 <- c(1, 1, 1, 1, 1, 1)
            sparse_bin_masses2 <- c(1, 1, 1, 1, 1, 1)
            sparse_bin_centres1 <- c(1, 2, 4, 7, 11, 16)
            sparse_bin_centres2 <- c(21, 22, 24, 27, 31, 36)
            sparse_histogram1 <- dhist(masses = sparse_bin_masses1, 
                                       locations = sparse_bin_centres1)
            sparse_histogram2 <- dhist(masses = sparse_bin_masses2, 
                                       locations = sparse_bin_centres2)
            
            dense_bin_centres1 <- 1:36
            dense_bin_centres2 <- dense_bin_centres1
            bin_mass_sequence <- c(1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1)
            bin_mass_padding <- rep(0, length(dense_bin_centres1) - 
                                      length(bin_mass_sequence))
            dense_bin_masses1 <- c(bin_mass_sequence, bin_mass_padding)
            dense_bin_masses2 <- c(bin_mass_padding, bin_mass_sequence)
            dense_histogram1 <- dhist(masses = dense_bin_masses1, 
                                      locations = dense_bin_centres1)
            dense_histogram2 <- dhist(masses = dense_bin_masses2, 
                                      locations = dense_bin_centres2)
            
            expect_equal(emd_lp(dense_bin_masses1, dense_bin_masses2,
                                dense_bin_centres1, dense_bin_centres2),
                         emd_lp(sparse_bin_masses1, sparse_bin_masses2,
                                sparse_bin_centres1, sparse_bin_centres2))
            expect_equal(emd_cs(dense_histogram1, dense_histogram2),
                         emd_cs(sparse_histogram1,sparse_histogram2))
            expect_equal(emd(dense_histogram1, dense_histogram2),
                         emd(sparse_histogram1, sparse_histogram2))
            })

test_that("EMD methods return same result when order of densely specified bins 
          is changed", {
            bin_masses1 <- c(1, 1, 1, 1, 0, 0, 0, 0, 0)
            bin_masses2 <- c(0, 0, 0, 0, 0, 1, 1, 1, 1)
            bin_centres1 <- 1:length(bin_masses1)
            bin_centres2 <- 1:length(bin_masses2)
            histogram1 <- dhist(masses = bin_masses1, locations = bin_centres1)
            histogram2 <- dhist(masses = bin_masses2, locations = bin_centres2)
            
            permuted_indexes1 <- sample(1:length(bin_masses1), replace = FALSE)
            permuted_indexes2 <- sample(1:length(bin_masses2), replace = FALSE)
            
            permuted_bin_masses1 <- bin_masses1[permuted_indexes1]
            permuted_bin_centres1 <- bin_centres1[permuted_indexes1]
            permuted_bin_masses2 <- bin_masses2[permuted_indexes2]
            permuted_bin_centres2 <- bin_centres2[permuted_indexes2]
            permuted_histogram1 <- dhist(masses = permuted_bin_masses1, 
                                         locations = permuted_bin_centres1)
            permuted_histogram2 <- dhist(masses = permuted_bin_masses2, 
                                         locations = permuted_bin_centres2)
            
            expect_equal(emd_lp(bin_masses1, bin_masses2, bin_centres1, bin_centres2),
                         emd_lp(permuted_bin_masses1, permuted_bin_masses2,
                                permuted_bin_centres1, permuted_bin_centres2))
            expect_equal(emd_cs(histogram1, histogram2),
                         emd_cs(permuted_histogram1, permuted_histogram2))
            expect_equal(emd(histogram1, histogram2), 
                         emd(permuted_histogram1, permuted_histogram2))
            })

test_that("EMD methods return same result when order of sparsely specified bins 
          is changed", {
            bin_masses1 <- c(1, 1, 1, 1, 1, 1)
            bin_masses2 <- c(1, 1, 1, 1, 1, 1)
            bin_centres1 <- c(1, 2, 4, 8, 16, 32)
            bin_centres2 <- c(-32, -16, -8, -4, -2, -1)
            histogram1 <- dhist(masses = bin_masses1, locations = bin_centres1)
            histogram2 <- dhist(masses = bin_masses2, locations = bin_centres2)
            
            permuted_indexes1 <- sample(1:length(bin_masses1), replace = FALSE)
            permuted_indexes2 <- sample(1:length(bin_masses2), replace = FALSE)
            
            permuted_bin_masses1 <- bin_masses1[permuted_indexes1]
            permuted_bin_centres1 <- bin_centres1[permuted_indexes1]
            permuted_bin_masses2 <- bin_masses2[permuted_indexes2]
            permuted_bin_centres2 <- bin_centres2[permuted_indexes2]
            permuted_histogram1 <- dhist(masses = permuted_bin_masses1, 
                                         locations = permuted_bin_centres1)
            permuted_histogram2 <- dhist(masses = permuted_bin_masses2, 
                                         locations = permuted_bin_centres2)
            
            expect_equal(emd_lp(bin_masses1, bin_masses2, bin_centres1, bin_centres2),
                         emd_lp(permuted_bin_masses1, permuted_bin_masses2,
                                permuted_bin_centres1, permuted_bin_centres2))
            expect_equal(emd_cs(histogram1, histogram2),
                         emd_cs(permuted_histogram1, permuted_histogram2))
            expect_equal(emd(histogram1, histogram2), 
                         emd(permuted_histogram1, permuted_histogram2))
          })

context("EMD: Next step")
test_that("next_step gives correct result for simple x1, x2", {
  x1 <- c(-3000, -2000, -1000, 0, 1000, 2000, 3000, 4000)
  x2 <- c(-3100, -2100, -1100, 10, 100, 1100, 2100, 2100, 4100)
  
  expected_shift <- 10
  actual_shift <- shift_to_next_alignment(x1, x2)
  expect_equal(actual_shift, expected_shift)
})

test_that("next_step gives correct result for random x1, x2", {
  test_fn <- function(n) {
    # Define x1 as random vector where minimum spacing between elements is
    # x1_prec
    x1_min <- -1000
    x1_max <- 1000
    x1_prec <- 10
    x1 <- unique(sort(trunc(runif(27, x1_min, x1_max)/x1_prec) * x1_prec))
    # Initialise x2 to a copy of x1 with all elements shifted right by 40% of
    # the minimum spacing between elements
    std_shift <- 0.4 * x1_prec
    x2 <- x1 + std_shift
    # Adjust a random element in x2 to be right-shifted by 10% rather than 40%
    min_shift <- 0.1 * x1_prec
    x2_rand_ind <- trunc(runif(1, 1, length(x2) + 1))
    x2[x2_rand_ind] <- x2[x2_rand_ind] - std_shift + min_shift
    
    expected_shift <- min_shift
    actual_shift <- shift_to_next_alignment(x1, x2)
    expect_equal(actual_shift, expected_shift)
  }
  purrr::walk(1:100, test_fn)
})

context("EMD: MinEMD exhaustive")
test_that("min_emd_exhaustive returns 0 when comparing a 1D feature distribution to itself",{
  bin_masses1 <- c(0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0)
  bin_masses2 <- bin_masses1
  bin_centres1 <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)
  bin_centres2 <- bin_centres1
  dhist1 <- dhist(masses = bin_masses1, locations = bin_centres1)
  dhist2 <- dhist(masses = bin_masses2, locations = bin_centres2)
  
  expected <- list(min_emd = 0, min_offset = 0)
  actual <- min_emd_exhaustive(dhist1, dhist2)
  expect_equal(actual, expected)
})

test_that("min_emd_exhaustive correct comparing a non-offset 1D feature distribution to itself",{
  bin_masses1 <- c(0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0)
  bin_masses2 <- bin_masses1
  bin_centres1 <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)
  bin_centres2 <- bin_centres1
  dhist1 <- dhist(masses = bin_masses1, locations = bin_centres1)
  dhist2 <- dhist(masses = bin_masses2, locations = bin_centres2)
  
  expected <- list(min_emd = 0, min_offset = 0)
  actual <- min_emd_exhaustive(dhist1, dhist2)
  expect_equal(actual, expected)
})

test_that("min_emd_exhaustive correct when comparing an offset 1D feature distribution to itself",{
  offset = 10
  bin_masses1 <- c(0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0)
  bin_masses2 <- bin_masses1
  bin_centres1 <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)
  bin_centres2 <- bin_centres1 + offset
  dhist1 <- dhist(masses = bin_masses1, locations = bin_centres1)
  dhist2 <- dhist(masses = bin_masses2, locations = bin_centres2)
  
  expected <- list(min_emd = 0, min_offset = offset)
  actual <- min_emd_exhaustive(dhist1, dhist2)
  expect_equal(actual, expected)
})
