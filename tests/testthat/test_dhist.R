library("purrr")
context("dhist: Discrete histogram from observations")

test_that("discrete_hist generates correct discrete histograms for random integer observations", {
  # Method for generating random observations containing specific locations a 
  # specific number of times
  random_observations <- function(locations, counts) {
    # Construct vector containing each location replicated "count" times
    observations <- purrr::simplify(purrr::map2(locations, counts, rep))
    # Randomise the order of the observations
    sample(observations, size = length(observations), replace = FALSE)
  }
  
  set.seed(2684)
  num_tests <- 100
  
  run_test <- function() {
    # Set parameters for generation of random observation sets
    num_observations <- 100
    location_range <- -(num_observations*3):(num_observations*3)
    # Do not allow zero counts as these locations will not be present in the 
    # observations generated from the locations and counts
    count_range <- 1:10
    
    # Generate random observation sets
    locations <- sample(location_range, num_observations, replace = FALSE)
    counts <- sample(count_range, num_observations, replace = TRUE)
    
    # Construct vector containing each location replicated "count" times
    observations_orig <- purrr::simplify(purrr::map2(locations, counts, rep))
    # Randomise the order of the observations
    observations <- sample(observations_orig, size = length(observations_orig), replace = FALSE)
    
    # Generate discrete histograms
    hist <- dhist_from_obs(observations)
    
    # discrete_hist will drop bins with zero counts, so remove these from the 
    # expected data (not necessary now we've restricted counts to be >= 1, but 
    # the bug where we generated test locations with zero counts was so annoying
    # to identify that we're going with a belt and braces approach)
    non_zero_count_indexes <- counts != 0
    expected_locations <- locations[non_zero_count_indexes]
    expected_counts <- counts[non_zero_count_indexes]
    # dhist_from_obs will return results with bins ordered by ascending location, 
    # so sort expected data to match
    sorted_locations <- sort(expected_locations, index.return = TRUE)
    sorted_location_indexes <- sorted_locations$ix
    expected_locations <- expected_locations[sorted_location_indexes]
    expected_counts <- expected_counts[sorted_location_indexes]
    
    # Check that histogram locations and counts match those used to generate the 
    # observations
    expect_true(all.equal(hist$locations, expected_locations))
    expect_true(all.equal(hist$masses, expected_counts))
  }
  
  for(i in 1:num_tests) {
    run_test()
  }
})

context("dhist: Discrete histogram variance")
test_that("dhist_variance returns sigma^2 for normal histograms", {
  num_hists <- 5
  num_bins <- 100001
  
  mus <- runif(num_hists, -10, 10)
  sigmas <- runif(num_hists, 0, 10)
  
  rand_locations <- function(mu, sigma) {return(seq(mu - 5 * sigma, mu + 5 * sigma, length.out = num_bins))}
  
  rand_dhists <- purrr::map2(mus, sigmas, function(mu, sigma) {
    locations <- rand_locations(mu, sigma)
    masses <- dnorm(locations, mean = mu, sd = sigma)
    return(dhist(masses = masses, locations = locations))
  })
  
  actuals <- purrr::map_dbl(rand_dhists, dhist_variance)
  expected <- purrr::map_dbl(sigmas, function(sigma) {return(sigma^2)})
  
  expect_equalish <- function(actual, expected) {
    scaled_diff <- abs(actual - expected)/min(actual, expected)
    max_diff <- 1e-4
    return(expect_lte(scaled_diff, max_diff))
  }
  purrr::map2(actuals, expected, expect_equalish)
})

context("dhist: Discrete histogram mass normalisation")
test_that("normalise_dhist_mass output sums to 1", {
  # Generate histograms with random masses (no centres needed for this test)
  num_hists <- 10
  num_bins <- 100
  
  mass_min <- 0
  mass_max <- 100
  rand_bin_masses <- function() {return(runif(num_bins, mass_min, mass_max))}
  bin_mass_lists <- replicate(num_hists, rand_bin_masses(), simplify = FALSE)

  actuals <- purrr::map(bin_mass_lists, function(masses) {
    # Locations are unimportant as they do not affect mass normalisation
    locations <- 1:length(masses)
    mass_normalised_dhist = normalise_dhist_mass(dhist(masses = masses, locations = locations))
    return(sum(mass_normalised_dhist$masses))
  })
  expected <- 1
  purrr::map_dbl(actuals, function(actual) {expect_equal(actual, expected)})
})

context("dhist: Discrete histogram variance normalisation")
test_that("normalise_histogram_variance output has variance of 1 for random histograms", {
  # Generate histograms with random masses and random centres
  num_hists <- 10
  num_bins <- 100
  
  mass_min <- 0
  mass_max <- 100
  rand_masses <- function() {return(runif(num_bins, mass_min, mass_max))}
  
  centre_min <- -30
  centre_max <- 70
  rand_locations <- function() {return(runif(num_bins, centre_min, centre_max))}
  
  rand_dhists <- replicate(num_hists, dhist(masses = rand_masses(), locations = rand_locations()), simplify = FALSE)
  
  actuals <- purrr::map(rand_dhists, function(dhist) {dhist_variance(normalise_dhist_variance(dhist))})
  expected <- 1
  purrr::map_dbl(actuals, function(actual) {expect_equal(actual, expected)})
})

test_that("normalise_histogram_variance output has variance of 1 for normal histograms", {
  num_hists <- 5
  num_bins <- 100001
  
  mus <- runif(num_hists, -10, 10)
  sigmas <- runif(num_hists, 0, 10)
  
  rand_locations <- function(mu, sigma) {return(seq(mu - 5 * sigma, mu + 5 * sigma, length.out = num_bins))}
  
  rand_dhists <- purrr::map2(mus, sigmas, function(mu, sigma) {
    locations <- rand_locations(mu, sigma)
    masses <- dnorm(locations, mean = mu, sd = sigma)
    return(dhist(masses = masses, locations = locations))
  })

  actuals <- purrr::map(rand_dhists, function(dhist) {dhist_variance(normalise_dhist_variance(dhist))})
  expected <- 1
  purrr::map_dbl(actuals, function(actual) {expect_equal(actual, expected)})
})

context("dhist: Sort dhist")
test_that("sort_dhist works", {
  # NOTE: Need to construct dhist objects explicitly as the dhist constructor
  # now returns a sorted dhist and we want to be independent of this
  dhist1 <- list(locations = c(7, 42, 1, 21, 101, 9), masses = c(15, 12, 16, 13, 11, 14))
  class(dhist1) <- "dhist"
  dhist2 <- list(locations = c(3, 0, -62, 7, 16, -58), masses = c(23, 24, 26, 22, 21, 25))
  class(dhist2) <- "dhist"
  
  expected1 = list(locations = c(1, 7, 9, 21, 42, 101), masses = c(16, 15, 14, 13 ,12, 11))
  class(expected1) <- "dhist"
  expected2 = list(locations = c(-62, -58, 0, 3, 7, 16), masses = c(26, 25, 24, 23, 22, 21))
  class(expected2) <- "dhist"
  
  actual1 <- sort_dhist(dhist1)
  actual2 <- sort_dhist(dhist2)

  expect_equal(actual1, expected1)
  expect_equal(actual2, expected2)
})

context("dhist: ECMF (norm_mass = F, norm_var = F)")
test_that("dhist_ecmf returns correct step function when smoothing_window_width is zero, normalise_mass is FALSE and normalise_variance is FALSE", {
  dhist1 <- dhist(locations = c(1, 2, 4, 7, 11, 16, 22), masses = c(21, 22, 23, 27, 31, 36, 42))
  
  ecmf1 <- dhist_ecmf(dhist1, smoothing_window_width = 0, normalise_mass = FALSE, normalise_variance = FALSE)
  actual_knots1 <- knots(ecmf1)
  actual_knots_ecds1 <- ecmf1(actual_knots1)
  inter_knots_x <- head(actual_knots1, length(actual_knots1) - 1)
  actual_inter_knots_ecds1 <- ecmf1(inter_knots_x)
  extra_knots <- c(actual_knots1[1] - 1, actual_knots1[length(actual_knots1)] + 1)
  actual_extra_knots_ecds1 <- ecmf1(extra_knots)
  
  cum_masses1 <- cumsum(dhist1$masses)
  max_cum_mass <- cum_masses1[length(cum_masses1)]
  expected_knots_ecds1 <- cum_masses1
  expected_inter_knots_ecds1 <- head(expected_knots_ecds1, length(expected_knots_ecds1) -1)
  expected_extra_knots_ecds1 <- c(0, max_cum_mass)
  expected_knots1 <- dhist1$locations
  
  expect_equal(actual_knots1, expected_knots1)
  expect_equal(actual_knots_ecds1, expected_knots_ecds1)
  expect_equal(actual_inter_knots_ecds1, expected_inter_knots_ecds1)
  expect_equal(actual_extra_knots_ecds1, expected_extra_knots_ecds1)
})
context("dhist: ECMF (norm_mass = T, norm_var = F)")
test_that("dhist_ecmf returns correct step function when smoothing_window_width is zero, normalise_mass is FALSE and normalise_variance is TRUE", {
  dhist1 <- dhist(locations = c(1, 2, 4, 7, 11, 16, 22), masses = c(21, 22, 23, 27, 31, 36, 42))
  
  ecmf1 <- dhist_ecmf(dhist1, smoothing_window_width = 0, normalise_mass = TRUE, normalise_variance = FALSE)
  actual_knots1 <- knots(ecmf1)
  actual_knots_ecds1 <- ecmf1(actual_knots1)
  inter_knots_x <- head(actual_knots1, length(actual_knots1) - 1)
  actual_inter_knots_ecds1 <- ecmf1(inter_knots_x)
  extra_knots <- c(actual_knots1[1] - 1, actual_knots1[length(actual_knots1)] + 1)
  actual_extra_knots_ecds1 <- ecmf1(extra_knots)
  
  normalised_cum_masses1 <- cumsum(dhist1$masses / sum(dhist1$masses))
  expected_knots_ecds1 <- normalised_cum_masses1
  expected_inter_knots_ecds1 <- head(expected_knots_ecds1, length(expected_knots_ecds1) -1)
  expected_extra_knots_ecds1 <- c(0,1)
  expected_knots1 <- dhist1$locations
  
  expect_equal(actual_knots1, expected_knots1)
  expect_equal(actual_knots_ecds1, expected_knots_ecds1)
  expect_equal(actual_inter_knots_ecds1, expected_inter_knots_ecds1)
  expect_equal(actual_extra_knots_ecds1, expected_extra_knots_ecds1)
})
context("dhist: ECMF (norm_mass = F, norm_var = T)")
test_that("dhist_ecdf returns correct step function when smoothing_window_width is zero, normalise_mass is FALSE and normalise_variance is TRUE", {
  dhist1 <- dhist(locations = c(1, 2, 4, 7, 11, 16, 22), masses = c(21, 22, 23, 27, 31, 36, 42))
  
  ecmf1 <- dhist_ecmf(dhist1, smoothing_window_width = 0, normalise_mass = FALSE, normalise_variance = TRUE)
  actual_knots1 <- knots(ecmf1)
  actual_knots_ecds1 <- ecmf1(actual_knots1)
  inter_knots_x <- head(actual_knots1, length(actual_knots1) - 1)
  actual_inter_knots_ecds1 <- ecmf1(inter_knots_x)
  extra_knots <- c(actual_knots1[1] - 1, actual_knots1[length(actual_knots1)] + 1)
  actual_extra_knots_ecds1 <- ecmf1(extra_knots)
  
  cum_masses1 <- cumsum(dhist1$masses)
  max_cum_mass <- cum_masses1[length(cum_masses1)]
  expected_knots_ecds1 <- cum_masses1
  expected_inter_knots_ecds1 <- head(expected_knots_ecds1, length(expected_knots_ecds1) -1)
  expected_extra_knots_ecds1 <- c(0, max_cum_mass)
  expected_knots1 <- normalise_dhist_variance(dhist1)$locations
  
  expect_equal(actual_knots1, expected_knots1)
  expect_equal(actual_knots_ecds1, expected_knots_ecds1)
  expect_equal(actual_inter_knots_ecds1, expected_inter_knots_ecds1)
  expect_equal(actual_extra_knots_ecds1, expected_extra_knots_ecds1)
})
context("dhist: ECMF (norm_mass = F, norm_var = T)")
test_that("dhist_ecdf returns correct step function when smoothing_window_width is zero and normalise_variance is true", {
  dhist1 <- dhist(locations = c(1, 2, 4, 7, 11, 16, 22), masses = c(21, 22, 23, 27, 31, 36, 42))
  
  ecmf1 <- dhist_ecmf(dhist1, smoothing_window_width = 0, normalise_mass = TRUE, normalise_variance = TRUE)
  actual_knots1 <- knots(ecmf1)
  actual_knots_ecds1 <- ecmf1(actual_knots1)
  inter_knots_x <- head(actual_knots1, length(actual_knots1) - 1)
  actual_inter_knots_ecds1 <- ecmf1(inter_knots_x)
  extra_knots <- c(actual_knots1[1] - 1, actual_knots1[length(actual_knots1)] + 1)
  actual_extra_knots_ecds1 <- ecmf1(extra_knots)
  
  normalised_cum_masses1 <- cumsum(dhist1$masses / sum(dhist1$masses))
  expected_knots_ecds1 <- normalised_cum_masses1
  expected_inter_knots_ecds1 <- head(expected_knots_ecds1, length(expected_knots_ecds1) -1)
  expected_extra_knots_ecds1 <- c(0,1)
  expected_knots1 <- normalise_dhist_variance(dhist1)$locations
  
  expect_equal(actual_knots1, expected_knots1)
  expect_equal(actual_knots_ecds1, expected_knots_ecds1)
  expect_equal(actual_inter_knots_ecds1, expected_inter_knots_ecds1)
  expect_equal(actual_extra_knots_ecds1, expected_extra_knots_ecds1)
})

context("dhist: Area between ECMFs")
test_that("area_between_dhist_ecmfs returns correct value", {
  # Example dhists constructed by hand to result in lots of "bowtie" segments   
  # for smoothed ECMFs and to allow expected areas to be calculated by hand
  dhistA <- dhist(locations = c(1, 3, 4), masses = c(2, 1, 1))
  dhistB <- dhist(locations = c(0, 2, 4, 5), masses = c(0.5, 2, 0.5, 1))
  expected_area_unsmoothed <- 4
  expected_area_smoothed <- 3
  
  # Generate ecmfs
  ecmfA_unsmoothed <- dhist_ecmf(dhistA, smoothing_window_width = 0)
  ecmfB_unsmoothed <- dhist_ecmf(dhistB, smoothing_window_width = 0)
  ecmfA_smoothed <- dhist_ecmf(dhistA, smoothing_window_width = 1)
  ecmfB_smoothed <- dhist_ecmf(dhistB, smoothing_window_width = 1)

  # Calculate area between ECMFs
  actual_area_unsmoothed <- area_between_dhist_ecmfs(ecmfA_unsmoothed, ecmfB_unsmoothed)
  actual_area_smoothed <- area_between_dhist_ecmfs(ecmfA_smoothed, ecmfB_smoothed)
  
  # Compare caculated areas with expected areas
  expect_equal(actual_area_unsmoothed, expected_area_unsmoothed)
  expect_equal(actual_area_smoothed, expected_area_smoothed)
})

context("dhist: Harmonise dhist locations")
test_that("harmonise_dhist_locations works A", {
  dhist1 <- dhist(masses = c(11, 12, 13), locations = c(1, 3, 5))
  dhist2 <- dhist(masses = c(21, 22, 23), locations = c(2, 4, 6))
  
  expected <- list(
    dhist1 = dhist(masses = c(11, 0, 12, 0, 13, 0), locations = c(1, 2, 3, 4, 5, 6)),
    dhist2 = dhist(masses = c(0, 21, 0, 22, 0, 23), locations = c(1, 2, 3, 4, 5, 6))
  )
  actual <- harmonise_dhist_locations(dhist1, dhist2)
  expect_equal(actual, expected)
})

test_that("harmonise_dhist_locations works B", {
  dhist1 <- dhist(masses = c(1, 1, 1), locations = c(1, 3, 5))
  dhist2 <- dhist(masses = c(1, 1, 1), locations = c(4, 5, 6))
  
  expected <- list(
    dhist1 = dhist(masses = c(1, 1, 1, 0, 0), locations = c(1, 3, 5, 4, 6)),
    dhist2 = dhist(masses = c(1, 1, 1, 0, 0), locations = c(4, 5, 6, 1, 3))
  )

  actual <- harmonise_dhist_locations(dhist1, dhist2)
  expect_equal(actual, expected)
})
