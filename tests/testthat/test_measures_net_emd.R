library("netdist")
library("purrr")
context("Measures NetEMD")


# BIN NORMALISATION: Property based tests
test_that("normalise_histogram_mass output sums to 1", {
  # Generate histograms with random masses (no centres needed for this test)
  num_hists <- 10
  num_bins <- 100
  
  mass_min <- 0
  mass_max <- 100
  rand_bin_masses <- function() {return(runif(num_bins, mass_min, mass_max))}
  bin_mass_lists <- replicate(num_hists, rand_bin_masses(), simplify = FALSE)
  
  actuals <- purrr::map_dbl(bin_mass_lists, function(bin_masses) {sum(normalise_histogram_mass(bin_masses))})
  expected <- 1
  purrr::map_dbl(actuals, function(actual) {expect_equal(actual, expected)})
})

test_that("histogram_variance returns sigma^2 for normal histograms", {
  num_hists <- 5
  num_bins <- 100001
  
  mus <- runif(num_hists, -10, 10)
  sigmas <- runif(num_hists, 0, 10)
  
  rand_bin_centres <- function(mu, sigma) {return(seq(mu - 5 * sigma, mu + 5 * sigma, length.out = num_bins))}
  
  bin_centre_lists <- purrr::map2(mus, sigmas, rand_bin_centres)
  bin_mass_lists <- purrr::pmap(list(bin_centre_lists, mus, sigmas), dnorm)
  
  actuals <- purrr::map2(bin_mass_lists, bin_centre_lists, histogram_variance)
  expected <- purrr::map(sigmas, function(sigma) {return(sigma^2)})
  
  expect_equalish <- function(actual, expected) {
    scaled_diff <- abs(actual - expected)/min(actual, expected)
    max_diff <- 1e-4
    return(expect_lte(scaled_diff, max_diff))
  }
  purrr::map2(actuals, expected, expect_equalish)
})

test_that("normalise_histogram_variance output has variance of 1 for normal histograms", {
  num_hists <- 5
  num_bins <- 100001
  
  mus <- runif(num_hists, -10, 10)
  sigmas <- runif(num_hists, 0, 10)
  
  rand_bin_centres <- function(mu, sigma) {return(seq(mu - 5 * sigma, mu + 5 * sigma, length.out = num_bins))}
  
  bin_centre_lists <- purrr::map2(mus, sigmas, rand_bin_centres)
  bin_mass_lists <- purrr::pmap(list(bin_centre_lists, mus, sigmas), dnorm)
  
  normalised_histogram_variance <- function(bin_masses, bin_centres) {
    histogram_variance(bin_masses, normalise_histogram_variance(bin_masses, bin_centres))
  }
  
  actuals <- purrr::map2_dbl(bin_mass_lists, bin_centre_lists, normalised_histogram_variance)
  expected <- 1
  purrr::map_dbl(actuals, function(actual) {expect_equal(actual, expected)})
})

test_that("normalise_histogram_variance output has variance of 1", {
  # Generate histograms with random masses and random centres
  num_hists <- 10
  num_bins <- 100
  
  mass_min <- 0
  mass_max <- 100
  rand_bin_masses <- function() {return(runif(num_bins, mass_min, mass_max))}
  bin_mass_lists <- replicate(num_hists, rand_bin_masses(), simplify = FALSE)
  
  centre_min <- -30
  centre_max <- 70
  rand_bin_centres <- function() {return(runif(num_bins, centre_min, centre_max))}
  bin_centre_lists <- replicate(num_hists, rand_bin_centres(), simplify = FALSE)
  
  normalised_histogram_variance <- function(bin_masses, bin_centres) {
    histogram_variance(bin_masses, normalise_histogram_variance(bin_masses, bin_centres))
  }
  
  actuals <- purrr::map2_dbl(bin_mass_lists, bin_centre_lists, normalised_histogram_variance)
  expected <- 1
  purrr::map_dbl(actuals, function(actual) {expect_equal(actual, expected)})
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

# EMD: Property-based tests
test_that("EMD methods return 0 when comparing a 1D feature distribution to itself",{
  bin_masses1 <- c(0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0)
  bin_masses2 <- bin_masses1
  bin_centres1 <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)
  bin_centres2 <- bin_centres1
  histogram1 <- list(masses = bin_masses1, locations = bin_centres1)
  histogram2 <- list(masses = bin_masses2, locations = bin_centres2)
  
  expected <- 0
  expect_equal(emd_lp(bin_masses1, bin_masses2, bin_centres1, bin_centres2), expected)
  expect_equal(emd_cs(bin_masses1, bin_masses2, bin_centres1, bin_centres2), expected)
  expect_equal(emd(histogram1, histogram2), expected)
})

test_that("EMD methods return numBins/2 when offsetting a symmetric discrete triangle distribution by 1", {
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
  histogram1 <- list(masses = bin_masses1, locations = bin_centres1)
  histogram2 <- list(masses = bin_masses2, locations = bin_centres2)
  num_nonzero_bins <- sum(bin_masses1 > 0)
  expected <- cost_fn(num_nonzero_bins)
  emd_lp(bin_masses1, bin_masses2, bin_centres1, bin_centres2)
  expect_equal(emd_lp(bin_masses1, bin_masses2, bin_centres1, bin_centres2), expected)
  expect_equal(emd_cs(bin_masses1, bin_masses2, bin_centres1, bin_centres2), expected)
  expect_equal(emd(histogram1, histogram2), expected)
  
  # Triangle(4, even), shifting by changing centres
  bin_masses1 <- c(0, 1, 2, 3, 4, 4, 3, 2, 1, 0)
  bin_masses2 <- c(0, 1, 2, 3, 4, 4, 3, 2, 1, 0)
  bin_centres1 <- 1:length(bin_masses1)
  bin_centres2 <- 1:length(bin_masses2) + 1
  histogram1 <- list(masses = bin_masses1, locations = bin_centres1)
  histogram2 <- list(masses = bin_masses2, locations = bin_centres2)
  num_nonzero_bins <- sum(bin_masses1 > 0)
  expected <- cost_fn(num_nonzero_bins)
  expect_equal(emd_lp(bin_masses1, bin_masses2, bin_centres1, bin_centres2), expected)
  expect_equal(emd_cs(bin_masses1, bin_masses2, bin_centres1, bin_centres2), expected)
  expect_equal(emd(histogram1, histogram2), expected)
  
  # Triangle(5, odd), shifting by changing masses
  bin_masses1 <- c(0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0)
  bin_masses2 <- c(0, 0, 1, 2, 3, 4, 5, 4, 3, 2, 1)
  bin_centres1 <- 1:length(bin_masses1)
  bin_centres2 <- 1:length(bin_masses2)
  histogram1 <- list(masses = bin_masses1, locations = bin_centres1)
  histogram2 <- list(masses = bin_masses2, locations = bin_centres2)
  expected <- cost_fn(sum(bin_masses1 > 0))
  expect_equal(emd_lp(bin_masses1, bin_masses2, bin_centres1, bin_centres2), expected)
  expect_equal(emd_cs(bin_masses1, bin_masses2, bin_centres1, bin_centres2), expected)
  expect_equal(emd(histogram1, histogram2), expected)
  
  # Triangle(5, odd), shifting by changing masses
  bin_masses1 <- c(0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0)
  bin_masses2 <- c(0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0)
  bin_centres1 <- 1:length(bin_masses1)
  bin_centres2 <- 1:length(bin_masses2) + 1
  histogram1 <- list(masses = bin_masses1, locations = bin_centres1)
  histogram2 <- list(masses = bin_masses2, locations = bin_centres2)
  expected <- cost_fn(sum(bin_masses1 > 0))
  expect_equal(emd_lp(bin_masses1, bin_masses2, bin_centres1, bin_centres2), expected)
  expect_equal(emd_cs(bin_masses1, bin_masses2, bin_centres1, bin_centres2), expected)
  expect_equal(emd(histogram1, histogram2), expected)
  
  # Triangle(5, even), shifting by changing masses
  bin_masses1 <- c(0, 1, 2, 3, 4, 5, 5, 4, 3, 2, 1, 0)
  bin_masses2 <- c(0, 0, 1, 2, 3, 4, 5, 5, 4, 3, 2, 1)
  bin_centres1 <- 1:length(bin_masses1)
  bin_centres2 <- 1:length(bin_masses2)
  histogram1 <- list(masses = bin_masses1, locations = bin_centres1)
  histogram2 <- list(masses = bin_masses2, locations = bin_centres2)
  expected <- cost_fn(sum(bin_masses1 > 0))
  expect_equal(emd_lp(bin_masses1, bin_masses2, bin_centres1, bin_centres2), expected)
  expect_equal(emd_cs(bin_masses1, bin_masses2, bin_centres1, bin_centres2), expected)
  expect_equal(emd(histogram1, histogram2), expected)
  
  # Triangle(5, even), shifting by changing centres
  bin_masses1 <- c(0, 1, 2, 3, 4, 5, 5, 4, 3, 2, 1, 0)
  bin_masses2 <- c(0, 1, 2, 3, 4, 5, 5, 4, 3, 2, 1, 0)
  bin_centres1 <- 1:length(bin_masses1)
  bin_centres2 <- 1:length(bin_masses2) + 1
  histogram1 <- list(masses = bin_masses1, locations = bin_centres1)
  histogram2 <- list(masses = bin_masses2, locations = bin_centres2)
  expected <- cost_fn(sum(bin_masses1 > 0))
  expect_equal(emd_lp(bin_masses1, bin_masses2, bin_centres1, bin_centres2), expected)
  expect_equal(emd_cs(bin_masses1, bin_masses2, bin_centres1, bin_centres2), expected)
  expect_equal(emd(histogram1, histogram2), expected)
  
  # Triangle(6, odd), shifting by changing masses
  bin_masses1 <- c(0, 1, 2, 3, 4, 5, 6, 5, 4, 3, 2, 1, 0)
  bin_masses2 <- c(0, 0, 1, 2, 3, 4, 5, 6, 5, 4, 3, 2, 1)
  bin_centres1 <- 1:length(bin_masses1)
  bin_centres2 <- 1:length(bin_masses2)
  histogram1 <- list(masses = bin_masses1, locations = bin_centres1)
  histogram2 <- list(masses = bin_masses2, locations = bin_centres2)
  expected <- cost_fn(sum(bin_masses1 > 0))
  expect_equal(emd_lp(bin_masses1, bin_masses2, bin_centres1, bin_centres2), expected)
  expect_equal(emd_cs(bin_masses1, bin_masses2, bin_centres1, bin_centres2), expected)
  expect_equal(emd(histogram1, histogram2), expected)
  
  # Triangle(6, odd), shifting by changing centres
  bin_masses1 <- c(0, 1, 2, 3, 4, 5, 6, 5, 4, 3, 2, 1, 0)
  bin_masses2 <- c(0, 1, 2, 3, 4, 5, 6, 5, 4, 3, 2, 1, 0)
  bin_centres1 <- 1:length(bin_masses1)
  bin_centres2 <- 1:length(bin_masses2) + 1
  histogram1 <- list(masses = bin_masses1, locations = bin_centres1)
  histogram2 <- list(masses = bin_masses2, locations = bin_centres2)
  expected <- cost_fn(sum(bin_masses1 > 0))
  expect_equal(emd_lp(bin_masses1, bin_masses2, bin_centres1, bin_centres2), expected)
  expect_equal(emd_cs(bin_masses1, bin_masses2, bin_centres1, bin_centres2), expected)
  expect_equal(emd(histogram1, histogram2), expected)
  
})

test_that("EMD methods return same result for densely and sparsely specified bins", {
  sparse_bin_masses1 <- c(1, 1, 1, 1, 1, 1)
  sparse_bin_masses2 <- c(1, 1, 1, 1, 1, 1)
  sparse_bin_centres1 <- c(1, 2, 4, 7, 11, 16)
  sparse_bin_centres2 <- c(21, 22, 24, 27, 31, 36)
  sparse_histogram1 <- list(masses = sparse_bin_masses1, locations = sparse_bin_centres1)
  sparse_histogram2 <- list(masses = sparse_bin_masses2, locations = sparse_bin_centres2)
  
  dense_bin_centres1 <- 1:36
  dense_bin_centres2 <- dense_bin_centres1
  bin_mass_sequence <- c(1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1)
  bin_mass_padding <- rep(0, length(dense_bin_centres1) - length(bin_mass_sequence))
  dense_bin_masses1 <- c(bin_mass_sequence, bin_mass_padding)
  dense_bin_masses2 <- c(bin_mass_padding, bin_mass_sequence)
  dense_histogram1 <- list(masses = dense_bin_masses1, locations = dense_bin_centres1)
  dense_histogram2 <- list(masses = dense_bin_masses2, locations = dense_bin_centres2)
  
  expect_equal(emd_lp(dense_bin_masses1, dense_bin_masses2, 
                    dense_bin_centres1, dense_bin_centres2),
              emd_lp(sparse_bin_masses1, sparse_bin_masses2, 
                    sparse_bin_centres1, sparse_bin_centres2))
  expect_equal(emd_cs(dense_bin_masses1, dense_bin_masses2, 
                    dense_bin_centres1, dense_bin_centres2),
              emd_cs(sparse_bin_masses1, sparse_bin_masses2, 
                    sparse_bin_centres1, sparse_bin_centres2))
  expect_equal(emd(dense_histogram1, dense_histogram2), 
               emd(sparse_histogram1, sparse_histogram2))
})

test_that("EMD methods return same result when order of densely specified bins is changed", {
  bin_masses1 <- c(1, 1, 1, 1, 0, 0, 0, 0, 0)
  bin_masses2 <- c(0, 0, 0, 0, 0, 1, 1, 1, 1)
  bin_centres1 <- 1:length(bin_masses1)
  bin_centres2 <- 1:length(bin_masses2)
  histogram1 <- list(masses = bin_masses1, locations = bin_centres1)
  histogram2 <- list(masses = bin_masses2, locations = bin_centres2)
  
  permuted_indexes1 <- sample(1:length(bin_masses1), replace = FALSE)
  permuted_indexes2 <- sample(1:length(bin_masses2), replace = FALSE)
  
  permuted_bin_masses1 <- bin_masses1[permuted_indexes1]
  permuted_bin_centres1 <- bin_centres1[permuted_indexes1]
  permuted_bin_masses2 <- bin_masses2[permuted_indexes2]
  permuted_bin_centres2 <- bin_centres2[permuted_indexes2]
  permuted_histogram1 <- list(masses = permuted_bin_masses1, locations = permuted_bin_centres1)
  permuted_histogram2 <- list(masses = permuted_bin_masses2, locations = permuted_bin_centres2)
  
  expect_equal(emd_lp(bin_masses1, bin_masses2, bin_centres1, bin_centres2),
              emd_lp(permuted_bin_masses1, permuted_bin_masses2, 
                    permuted_bin_centres1, permuted_bin_centres2))
  expect_equal(emd_cs(bin_masses1, bin_masses2, bin_centres1, bin_centres2),
              emd_cs(permuted_bin_masses1, permuted_bin_masses2, 
                    permuted_bin_centres1, permuted_bin_centres2))
  expect_equal(emd(histogram1, histogram2), emd(permuted_histogram1, permuted_histogram2))
})

test_that("EMD methods return same result when order of sparsely specified bins is changed", {
  bin_masses1 <- c(1, 1, 1, 1, 1, 1)
  bin_masses2 <- c(1, 1, 1, 1, 1, 1)
  bin_centres1 <- c(1, 2, 4, 8, 16, 32)
  bin_centres2 <- c(-32, -16, -8, -4, -2, -1)
  histogram1 <- list(masses = bin_masses1, locations = bin_centres1)
  histogram2 <- list(masses = bin_masses2, locations = bin_centres2)
  
  permuted_indexes1 <- sample(1:length(bin_masses1), replace = FALSE)
  permuted_indexes2 <- sample(1:length(bin_masses2), replace = FALSE)
  
  permuted_bin_masses1 <- bin_masses1[permuted_indexes1]
  permuted_bin_centres1 <- bin_centres1[permuted_indexes1]
  permuted_bin_masses2 <- bin_masses2[permuted_indexes2]
  permuted_bin_centres2 <- bin_centres2[permuted_indexes2]
  permuted_histogram1 <- list(masses = permuted_bin_masses1, locations = permuted_bin_centres1)
  permuted_histogram2 <- list(masses = permuted_bin_masses2, locations = permuted_bin_centres2)
  
  expect_equal(emd_lp(bin_masses1, bin_masses2, bin_centres1, bin_centres2),
               emd_lp(permuted_bin_masses1, permuted_bin_masses2, 
                      permuted_bin_centres1, permuted_bin_centres2))
  expect_equal(emd_cs(bin_masses1, bin_masses2, bin_centres1, bin_centres2),
               emd_cs(permuted_bin_masses1, permuted_bin_masses2, 
                      permuted_bin_centres1, permuted_bin_centres2))
  expect_equal(emd(histogram1, histogram2), emd(permuted_histogram1, permuted_histogram2))
  })

# NetEMD: Property-based tests
test_that("net_emd returns 0 when comparing an integer location histogram against itself", {
  
  self_net_emd <- function(histogram, shift) {
    net_emd(histogram, shift_dhist(histogram, shift))
  }
  expected <- 0
  
  locations <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)
  masses <- c(0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0)
  histogram <- dhist(locations = locations, masses = masses)

  expect_equal(self_net_emd(histogram, shift = 1), expected)
  expect_equal(self_net_emd(histogram, shift = 0.5), expected)
  expect_equal(self_net_emd(histogram, shift = 0.1), expected)
  expect_equal(self_net_emd(histogram, shift = 0.05), expected)
  expect_equal(self_net_emd(histogram, shift = 0.01), expected)
})

test_that("net_emd returns 0 when comparing any normal histogram against itself (no offset)", {
  num_hists <- 5
  num_bins <- 101
  
  mus <- runif(num_hists, -10, 10)
  sigmas <- runif(num_hists, 0, 10)

  rand_locations <- function(mu, sigma) {return(seq(mu - 5 * sigma, mu + 5 * sigma, length.out = num_bins))}
  
  location_lists <- purrr::map2(mus, sigmas, rand_locations)
  mass_lists <- purrr::pmap(list(location_lists, mus, sigmas), dnorm)
  
  self_net_emd <- function(locations, masses) {
    histogram <- dhist(locations = locations, masses = masses)
    net_emd(histogram, histogram)
  }
  
  expect_equalish <- function(actual, expected) {
    diff <- abs(actual - expected)
    max_diff <- 1e-12
    return(diff <= max_diff)
  }
  
  expected <- 0
  actuals <- purrr::map2_dbl(location_lists, mass_lists, self_net_emd)
  purrr::map_dbl(actuals, function(actual) {expect_equalish(actual, expected)})

})

test_that("net_emd returns 0 when comparing any normal histogram randomly offset against itself", {
  
  num_hists <- 2
  num_bins <- 101
  num_offsets <- 3

  mus <- runif(num_hists, -10, 10)
  sigmas <- runif(num_hists, 0, 10)
  offsets <- runif(num_offsets, -10, 10)

  rand_locations <- function(mu, sigma) {return(seq(mu - 5 * sigma, mu + 5 * sigma, length.out = num_bins))}
  
  location_lists <- purrr::map2(mus, sigmas, rand_locations)
  mass_lists <- purrr::pmap(list(location_lists, mus, sigmas), dnorm)
  offset_lists <- replicate(num_hists, offsets, simplify = FALSE)
  
  expect_equalish <- function(actual, expected) {
    diff <- abs(actual - expected)
    max_diff <- 1e-12
    return(diff <= max_diff)
  }
  
  net_emd_offset_self <- function(locations, masses, offsets) {
    histogram <- dhist(locations = locations, masses = masses)
    net_emds <- purrr::map_dbl(offsets, function(offset) {net_emd(histogram, shift_dhist(histogram, offset))})
    return(net_emds)
  }

  expected <- 0
  actuals_list <- purrr::pmap(list(location_lists, mass_lists, offset_lists), net_emd_offset_self)
  purrr::map(actuals_list, function(actuals) {
        purrr::map_dbl(actuals, function(actual) {expect_equalish(actual, expected)})
  })
})

# EMD and NET_EMD: Virus PPI datasets
test_that("emd return 0 when comparing graphlet orbit degree distributions of virus PPI graphs to themselves", {
  # Load viurs PPI network data in ORCA-compatible edge list format
  data("virusppi")
  data_indexes <- 1:length(virusppi)
  data_names <- attr(virusppi, "name")

  # Calculate graphlet orbit degree distributions up to graphlet order 4
  virus_godd <- purrr::map(virusppi, godd)

  # Map over virus PPI networks
  purrr::walk(virus_godd, function(godd) {
    purrr::walk(godd, function(godd_Ox) {
      expect_equal(emd(godd_Ox, godd_Ox), 0)
    })
  })
})

test_that("net_emd return 0 when comparing graphlet orbit degree distributions of virus PPI graphs to themselves", {
  # Load viurs PPI network data in ORCA-compatible edge list format
  data("virusppi")
  data_indexes <- 1:length(virusppi)
  data_names <- attr(virusppi, "name")
  
  # Calculate graphlet orbit degree distributions up to graphlet order 4
  virus_godd <- purrr::map(virusppi, godd)
  
  # Map over virus PPI networks
  purrr::walk(virus_godd, function(godd) {
    purrr::walk(godd, function(godd_Ox) {
      expect_equal(net_emd(godd_Ox, godd_Ox), 0)
    })
  })
})
