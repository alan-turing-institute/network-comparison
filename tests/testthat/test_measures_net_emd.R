context("Measures NetEMD: Cost matrix")
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

context("Measures NetEMD: EMD")
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

context("Measures NetEMD: NetEMD")
# NetEMD: Property-based tests
test_that("net_emd returns 0 when comparing an integer location histogram offset
          against itself", {

  self_net_emd <- function(histogram, shift, method, step_size = NULL) {
    missing(step_size)
    net_emd(histogram, shift_dhist(histogram, shift), method, step_size)
  }
  expected <- 0

  locations <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)
  masses <- c(0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0)
  histogram <- dhist(locations = locations, masses = masses)

  expect_equal(self_net_emd(histogram, shift = 1, "optimise"), expected)
  expect_equal(self_net_emd(histogram, shift = 1, "fixed_step"), expected)
  expect_equal(self_net_emd(histogram, shift = 1, "fixed_step", step_size = 1), 
               expected)
  expect_equal(self_net_emd(histogram, shift = 0.5, "optimise"), expected)
  expect_equal(self_net_emd(histogram, shift = 0.5, "fixed_step"), expected)
  expect_equal(self_net_emd(histogram, shift = 1, "fixed_step", step_size = 1), 
               expected)
  expect_equal(self_net_emd(histogram, shift = 0.1, "optimise"), expected)
  expect_equal(self_net_emd(histogram, shift = 0.1, "fixed_step"), expected)
  expect_equal(self_net_emd(histogram, shift = 1, "fixed_step", step_size = 1), 
               expected)
  expect_equal(self_net_emd(histogram, shift = 0.05, "optimise"), expected)
  expect_equal(self_net_emd(histogram, shift = 0.05, "fixed_step"), expected)
  expect_equal(self_net_emd(histogram, shift = 1, "fixed_step", step_size = 1), 
               expected)
  expect_equal(self_net_emd(histogram, shift = 0.01, "optimise"), expected)
  expect_equal(self_net_emd(histogram, shift = 0.01, "fixed_step"), expected)
  expect_equal(self_net_emd(histogram, shift = 1, "fixed_step", step_size = 1), 
               expected)
  expect_equal(self_net_emd(histogram, shift = 0, "optimise"), expected)
  expect_equal(self_net_emd(histogram, shift = 0, "fixed_step"), expected)
  expect_equal(self_net_emd(histogram, shift = 1, "fixed_step", step_size = 1), 
               expected)
})

test_that("net_emd returns min_emd = 0 and min_offset = 0 when comparing an 
          integer location histogram offset against itself", {

  expect_self_net_emd_correct <- function(histogram, shift, method, 
                                          step_size = NULL, 
                                          return_details = FALSE) {
    self_net_emd <- net_emd(histogram, shift_dhist(histogram, shift), 
                            method, step_size, return_details)
    expected <- list(net_emd = 0, min_emds = 0, min_offsets = shift)
    expect_equal(self_net_emd, expected)
  }

  locations <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)
  masses <- c(0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0)
  histogram <- dhist(locations = locations, masses = masses)

  expect_self_net_emd_correct(histogram, shift = 1, "optimise", 
                              return_details = TRUE)
  expect_self_net_emd_correct(histogram, shift = 1, "fixed_step", 
                              return_details = TRUE)
  expect_self_net_emd_correct(histogram, shift = 1, "fixed_step", 
                              step_size = 1, return_details = TRUE)
  expect_self_net_emd_correct(histogram, shift = 0.5, "optimise", 
                              return_details = TRUE)
  expect_self_net_emd_correct(histogram, shift = 0.5, "fixed_step", 
                              return_details = TRUE)
  expect_self_net_emd_correct(histogram, shift = 1, "fixed_step", 
                              step_size = 1, return_details = TRUE)
  expect_self_net_emd_correct(histogram, shift = 0.1, "optimise", 
                              return_details = TRUE)
  expect_self_net_emd_correct(histogram, shift = 0.1, "fixed_step", 
                              return_details = TRUE)
  expect_self_net_emd_correct(histogram, shift = 1, "fixed_step", 
                              step_size = 1, return_details = TRUE)
  expect_self_net_emd_correct(histogram, shift = 0.05, "optimise", 
                              return_details = TRUE)
  expect_self_net_emd_correct(histogram, shift = 0.05, "fixed_step", 
                              return_details = TRUE)
  expect_self_net_emd_correct(histogram, shift = 1, "fixed_step", 
                              step_size = 1, return_details = TRUE)
  expect_self_net_emd_correct(histogram, shift = 0.01, "optimise", 
                              return_details = TRUE)
  expect_self_net_emd_correct(histogram, shift = 0.01, "fixed_step", 
                              return_details = TRUE)
  expect_self_net_emd_correct(histogram, shift = 1, "fixed_step", 
                              step_size = 1, return_details = TRUE)
  expect_self_net_emd_correct(histogram, shift = 0, "optimise", 
                              return_details = TRUE)
  expect_self_net_emd_correct(histogram, shift = 0, "fixed_step", 
                              return_details = TRUE)
  expect_self_net_emd_correct(histogram, shift = 1, "fixed_step", 
                              step_size = 1, return_details = TRUE)
})

test_that("net_emd returns 0 when comparing any normal histogram against itself (no offset)", {
  num_hists <- 5
  num_bins <- 101

  mus <- runif(num_hists, -10, 10)
  sigmas <- runif(num_hists, 0, 10)

  rand_locations <- function(mu, sigma) {
    return(seq(mu - 5 * sigma, mu + 5 * sigma, length.out = num_bins))
    }

  rand_dhists <- purrr::map2(mus, sigmas, function(mu, sigma) {
    locations <- rand_locations(mu, sigma)
    masses <- dnorm(locations, mean = mu, sd = sigma)
    return(dhist(masses = masses, locations = locations))
  })

  expected <- 0
  actuals_opt <- purrr::map(rand_dhists, function(dhist) {
    net_emd(dhist, dhist, method = "optimise")
    })
  purrr::map_dbl(actuals_opt, function(actual) {expect_equal(actual, expected)})

  actuals_step_default <- purrr::map(rand_dhists, function(dhist) {
    net_emd(dhist, dhist, method = "fixed_step")
    })
  purrr::map_dbl(actuals_step_default, function(actual) {
    expect_equal(actual, expected)
    })
})

test_that("net_emd returns 0 when comparing any normal histogram randomly offset
          against itself", {

  num_hists <- 2
  num_bins <- 101
  num_offsets <- 3

  mus <- runif(num_hists, -10, 10)
  sigmas <- runif(num_hists, 0, 10)
  offsets <- runif(num_offsets, -10, 10)

  rand_locations <- function(mu, sigma) {
    return(seq(mu - 5 * sigma, mu + 5 * sigma, length.out = num_bins))
    }

  rand_dhists <- purrr::map2(mus, sigmas, function(mu, sigma) {
    locations <- rand_locations(mu, sigma)
    masses <- dnorm(locations, mean = mu, sd = sigma)
    return(dhist(masses = masses, locations = locations))
  })

  offset_lists <- replicate(num_hists, offsets, simplify = FALSE)

  net_emd_offset_self <- function(dhist, offsets, method) {
    net_emds <- purrr::map_dbl(offsets, function(offset) {
      net_emd(dhist, shift_dhist(dhist, offset), method = method)})
    return(net_emds)
  }

  expected <- 0
  actuals_list_opt <- purrr::map2(rand_dhists, offset_lists, 
                                  function(dhist, offsets) {
    net_emd_offset_self(dhist, offsets, method = "optimise")})
  purrr::map(actuals_list_opt, function(actuals) {
        purrr::map_dbl(actuals, function(actual) {
          expect_equal(actual, expected)})
  })
  actuals_list_step <- purrr::map2(rand_dhists, offset_lists, 
                                   function(dhist, offsets) {
    net_emd_offset_self(dhist, offsets, method = "fixed_step")})
  purrr::map(actuals_list_step, function(actuals) {
    purrr::map_dbl(actuals, function(actual) {expect_equal(actual, expected)})
  })
})

test_that("net_emd returns min_emd = 0 and min_offset = 0 when comparing any 
          normal histogram randomly offset against itself", {

  num_hists <- 2
  num_bins <- 101
  num_offsets <- 3

  mus <- runif(num_hists, -10, 10)
  sigmas <- runif(num_hists, 0, 10)
  offsets <- runif(num_offsets, -10, 10)

  rand_locations <- function(mu, sigma) {
    return(seq(mu - 5 * sigma, mu + 5 * sigma, length.out = num_bins))}

  rand_dhists <- purrr::map2(mus, sigmas, function(mu, sigma) {
    locations <- rand_locations(mu, sigma)
    masses <- dnorm(locations, mean = mu, sd = sigma)
    return(dhist(masses = masses, locations = locations))
  })

  offset_lists <- replicate(num_hists, offsets, simplify = FALSE)

  expect_self_net_emd_correct <- 
    function(histogram, shift, method, step_size = NULL, 
             return_details = FALSE) {
    self_net_emd <- net_emd(histogram, shift_dhist(histogram, shift), 
                            method, step_size, return_details)
    expected <- list(net_emd = 0, min_emds = 0, min_offsets = shift)
    expect_equal(self_net_emd, expected)
  }

  purrr::map2(rand_dhists, offset_lists, function(dhist, offsets) {
    purrr::map(offsets, function(offset){
      expect_self_net_emd_correct(dhist, offset, "optimise", 
                                  return_details = TRUE)
    })
  })
  
  purrr::map2(rand_dhists, offset_lists, function(dhist, offsets) {
    purrr::map(offsets, function(offset){
      expect_self_net_emd_correct(dhist, offset, "fixed_step", r
                                  eturn_details = TRUE)
    })
  })
})

context("Measures NetEMD: Virus PPI (EMD)")
# EMD and NET_EMD: Virus PPI datasets
test_that("emd return 0 when comparing graphlet orbit degree distributions of 
          virus PPI graphs to themselves", {
  # Load viurs PPI network data in ORCA-compatible edge list format
  data("virusppi")
  data_indexes <- 1:length(virusppi)
  data_names <- attr(virusppi, "name")

  # Calculate graphlet-based degree distributions up to graphlet order 4
  virus_gdd <- purrr::map(virusppi, gdd)

  # Map over virus PPI networks
  purrr::walk(virus_gdd, function(gdd) {
    purrr::walk(gdd, function(gdd_Ox) {
      expect_equal(emd(gdd_Ox, gdd_Ox), 0)
    })
  })
})

context("Measures NetEMD: Virus PPI (NetEMD)")
test_that("net_emd return 0 when comparing graphlet orbit degree distributions 
          of virus PPI graphs to themselves", {
  # Load virus PPI network data in ORCA-compatible edge list format
  data("virusppi")
  data_indexes <- 1:length(virusppi)
  data_names <- attr(virusppi, "name")

  # Calculate graphlet-based degree distributions up to graphlet order 4
  virus_gdd <- purrr::map(virusppi, gdd)

  expect_equalish <- function(actual, expected) {
    diff <- abs(actual - expected)
    max_diff <- 1e-12
    return(expect_lte(diff, max_diff))
  }

  # Map over virus PPI networks
  purrr::walk(virus_gdd, function(gdd) {
    purrr::walk(gdd, function(gdd_Ox) {
      expect_equalish(net_emd(gdd_Ox, gdd_Ox, method = "optimise", 
                              smoothing_window_width = 0), 0)
    })
  })
})

context("Measures NetEMD: Random graphs (EMD)")
# EMD and NET_EMD: Random graph datasets
test_that("emd return 0 when comparing graphlet orbit degree distributions of 
          random graphs to themselves", {
  # Load random graph data in ORCA-compatible edge list format
  random_graphs <- read_all_graphs_as_orca_graphs(
    system.file(package = "netdist", "extdata", "random"),
    format = "ncol", pattern = "*")
  data_indexes <- 1:length(random_graphs)
  data_names <- attr(random_graphs, "name")

  # Calculate graphlet-based degree distributions up to graphlet order 4
  random_gdd <- purrr::map(random_graphs, gdd)

  # Map over random graphs
  purrr::walk(random_gdd, function(gdd) {
    purrr::walk(gdd, function(gdd_Ox) {
      expect_equal(emd(gdd_Ox, gdd_Ox), 0)
    })
  })
})

context("Measures NetEMD: Random graphs (NetEMD)")
test_that("net_emd return 0 when comparing graphlet orbit degree distributions 
          of random graphs to themselves", {
  # Load random graph data in ORCA-compatible edge list format
  random_graphs <- read_all_graphs_as_orca_graphs(
    system.file(package = "netdist", "extdata", "random"),
    format = "ncol", pattern = "*")
  data_indexes <- 1:length(random_graphs)
  data_names <- attr(random_graphs, "name")

  # Calculate graphlet-based degree distributions up to graphlet order 4
  random_gdd <- purrr::map(random_graphs, gdd)

  expect_equalish <- function(actual, expected) {
    diff <- abs(actual - expected)
    max_diff <- 1e-12
    return(expect_lte(diff, max_diff))
  }

  # Map over random graphs
  purrr::walk(random_gdd, function(gdd) {
    purrr::walk(gdd, function(gdd_Ox) {
      expect_equalish(net_emd(gdd_Ox, gdd_Ox, method = "optimise", 
                              smoothing_window_width = 0), 0)
    })
  })
})

context("Measures NetEMD: All graphs in directory")
test_that("net_emds_for_all_graphs works", {
  # Set source directory and file properties for Virus PPI graph edge files
  source_dir <- system.file(file.path("extdata", "VRPINS"), package = "netdist")
  edge_format = "ncol"
  file_pattern = ".txt"
  
  # Set number of threads to use at once for parallel processing.
  num_threads = getOption("mc.cores", 2L)
  
  # Use previously tested GDD code to generate inputs to function under test
  gdds_orbits_g4 <- gdd_for_all_graphs(
    source_dir = source_dir, format = edge_format, pattern = file_pattern, 
    feature_type = "orbit", max_graphlet_size = 4)
  gdds_orbits_g5 <- gdd_for_all_graphs(
    source_dir = source_dir, format = edge_format, pattern = file_pattern, 
    feature_type = "orbit", max_graphlet_size = 5)
  gdds_graphlets_g4 <- gdd_for_all_graphs(
    source_dir = source_dir, format = edge_format, pattern = file_pattern, 
    feature_type = "graphlet", max_graphlet_size = 4)
  gdds_graphlets_g5 <- gdd_for_all_graphs(
    source_dir = source_dir, format = edge_format, pattern = file_pattern, 
    feature_type = "graphlet", max_graphlet_size = 5)
  gdds_graphlets_g4_e1 <- gdd_for_all_graphs(
    source_dir = source_dir, format = edge_format, pattern = file_pattern, 
    feature_type = "graphlet", max_graphlet_size = 4, ego_neighbourhood_size = 1)
  gdds_graphlets_g5_e1 <- gdd_for_all_graphs(
    source_dir = source_dir, format = edge_format, pattern = file_pattern, 
    feature_type = "graphlet", max_graphlet_size = 5, ego_neighbourhood_size = 1)
  gdds_graphlets_g4_e2 <- gdd_for_all_graphs(
    source_dir = source_dir, format = edge_format, pattern = file_pattern, 
    feature_type = "graphlet", max_graphlet_size = 4, ego_neighbourhood_size = 2)
  gdds_graphlets_g5_e2 <- gdd_for_all_graphs(
    source_dir = source_dir, format = edge_format, pattern = file_pattern, 
    feature_type = "graphlet", max_graphlet_size = 5, ego_neighbourhood_size = 2)
  
  # Use previously tested NetEMD function to generate expected NetEMD scores
  # individually and combine into expected output for code under test
  expected_net_emd_fn<- function(gdds) {
    list(net_emds = c(net_emd(gdds$EBV, gdds$ECL), net_emd(gdds$EBV, gdds$HSV),
                      net_emd(gdds$EBV, gdds$KSHV), net_emd(gdds$EBV, gdds$VZV),
                      net_emd(gdds$ECL, gdds$HSV), net_emd(gdds$ECL, gdds$KSHV), 
                      net_emd(gdds$ECL, gdds$VZV), net_emd(gdds$HSV, gdds$KSHV), 
                      net_emd(gdds$HSV, gdds$VZV), net_emd(gdds$KSHV, gdds$VZV)),
         comp_spec = cross_comparison_spec(gdds))
  }
  
  # Comparison function for clarity
  compare_fn <- function(gdds) {
    expect_equal(net_emds_for_all_graphs(gdds), expected_net_emd_fn(gdds))
  }
  
  # Map over test parameters, comparing actual gdds to expected
  # No ego-networks
  compare_fn(gdds_orbits_g4)
  compare_fn(gdds_orbits_g5)
  compare_fn(gdds_graphlets_g4)
  compare_fn(gdds_graphlets_g5)
  # Ego networks of order 1
  compare_fn(gdds_graphlets_g4_e1)
  compare_fn(gdds_graphlets_g5_e1)
  # Ego networks of order 2
  compare_fn(gdds_graphlets_g4_e2)
  compare_fn(gdds_graphlets_g5_e2)
})

