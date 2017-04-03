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

context("dhist: constructor and as_* transformation functions")
test_that("dhist constuctor has correct locations and masses, sorted by ascending location", {
  locations1 = c(7, 42, 1, 21, 101, 9)
  masses1 = c(15, 12, 16, 13, 11, 14)
  actual1 <- dhist(locations = locations1, masses = masses1)
  locations2 = c(3, 0, -62, 7, 16, -58)
  masses2 = c(23, 24, 26, 22, 21, 25)
  actual2 <- dhist(locations = locations2, masses = masses2)
  
  expected_class <- "dhist"
  expected_smoothing_window_width <- 0
  
  expected1 = list(locations = c(1, 7, 9, 21, 42, 101), 
                   masses = c(16, 15, 14, 13 ,12, 11), 
                   smoothing_window_width = expected_smoothing_window_width)
  class(expected1) <- expected_class
  
  expected2 = list(locations = c(-62, -58, 0, 3, 7, 16), 
                   masses = c(26, 25, 24, 23, 22, 21), 
                   smoothing_window_width = expected_smoothing_window_width)
  class(expected2) <- expected_class
  
  expect_equal(actual1, expected1)
  expect_equal(actual2, expected2)
})

test_that("as_smoothed_dhist sets smoothing_window_width correctly", {
  dhist_pre <- dhist(locations <- c(7, 42, 1, 21, 101, 9), 
                  masses = c(15, 12, 16, 13, 11, 14))
  expected_smoothing_window_width_pre <- 0
  expected_smoothing_window_width_post <- 1
  
  expect_equal(dhist_pre$smoothing_window_width, 
               expected_smoothing_window_width_pre)
  dhist_post <- as_smoothed_dhist(dhist_pre, expected_smoothing_window_width_post)
  expect_equal(dhist_post$smoothing_window_width, 
               expected_smoothing_window_width_post)
})

test_that("as_unsmoothed_dhist sets smoothing_window_width correctly", {
  dhist_pre <- dhist(locations <- c(7, 42, 1, 21, 101, 9), 
                     masses = c(15, 12, 16, 13, 11, 14),
                     smoothing_window_width <- 1)
  expected_smoothing_window_width_pre <- 1
  expected_smoothing_window_width_post <- 0
  
  expect_equal(dhist_pre$smoothing_window_width, 
               expected_smoothing_window_width_pre)
  dhist_post <- as_smoothed_dhist(dhist_pre, expected_smoothing_window_width_post)
  expect_equal(dhist_post$smoothing_window_width, 
               expected_smoothing_window_width_post)
})

context("dhist: Discrete histogram variance")
test_that("dhist_variance difference for smoothed and unsmoothed dhists is smoothing_window_width^2 / 12", {
  dhist <- dhist(locations <- c(7, 42, 1, 21, 101, 9),  masses = c(15, 12, 16, 13, 11, 14))
  # Be careful: ensure that no smoothing window width results in overlapping bins
  smoothing_window_width_A <- 1
  smoothing_window_width_B <- 2
  dhist_unsmoothed <- as_unsmoothed_dhist(dhist)
  dhist_smoothed_A <- as_smoothed_dhist(dhist, smoothing_window_width_A)
  dhist_smoothed_B <- as_smoothed_dhist(dhist, smoothing_window_width_B)
  
  var_unsmoothed <- dhist_variance(dhist_unsmoothed)
  var_smoothed_A <- dhist_variance(dhist_smoothed_A)
  var_smoothed_B <- dhist_variance(dhist_smoothed_B)
  
  expected_var_smoothed_A <- var_unsmoothed + ((smoothing_window_width_A^2) / 12)
  expected_var_smoothed_B <- var_unsmoothed + ((smoothing_window_width_B^2) / 12)
  
  expect_equal(var_smoothed_A, expected_var_smoothed_A)
  expect_equal(var_smoothed_B, expected_var_smoothed_B)
})

test_that("dhist_variance returns sigma^2 for unsmoothed normal histograms", {
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
test_that("normalise_histogram_variance output has variance of 1 for random integer histograms", {
  # Generate histograms with random masses and random centres
  num_hists <- 10
  num_bins <- 70
  
  mass_min <- 0
  mass_max <- 100
  rand_masses <- function() {return(runif(num_bins, mass_min, mass_max))}
  
  centre_min <- -30
  centre_max <- 70
  rand_locations <- function() {return(round(sample(centre_min:centre_max, num_bins), digits = 0))}
  
  rand_dhists <- replicate(num_hists, dhist(masses = rand_masses(), locations = rand_locations()), simplify = FALSE)
  
  smoothing_window_width <- 1
  rand_dhists_unsmoothed <- purrr::map(rand_dhists, as_unsmoothed_dhist)
  rand_dhists_smoothed <- purrr::map(rand_dhists, as_smoothed_dhist, smoothing_window_width = smoothing_window_width)
  
  actual_unsmoothed <- purrr::map(rand_dhists_unsmoothed, function(dhist) {dhist_variance(normalise_dhist_variance(dhist))})
  actual_smoothed <- purrr::map(rand_dhists_smoothed, function(dhist) {dhist_variance(normalise_dhist_variance(dhist))})
  expected <- 1
  purrr::map_dbl(actual_unsmoothed, function(actual) {expect_equal(actual, expected)})
  purrr::map_dbl(actual_smoothed, function(actual) {expect_equal(actual, expected)})
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

context("dhist: ECMF")
test_that("dhist_ecmf returns correct step function when smoothing_window_width is zero", {
  dhist1 <- dhist(locations = c(1, 2, 4, 7, 11, 16, 22), masses = c(21, 22, 23, 27, 31, 36, 42))
  dhist1_unsmoothed <- as_unsmoothed_dhist(dhist1)
  
  ecmf1 <- dhist_ecmf(dhist1)
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

context("dhist: Area between ECMFs (simple integer dhists)")
test_that("area_between_dhist_ecmfs returns correct value for simple integer dhists", {
  # Example dhists constructed by hand to result in lots of "bowtie" segments   
  # for smoothed ECMFs and to allow expected areas to be calculated by hand
  # Unsmoothed locations are on an integer grid, smoothed bin edges are on a
  # half-integer grid
  # Smoothed and unsmoothed ECMF cumulative masses are on integer grid
  # Smoothed ECMF crossing points are on a quarter-integer grid
  dhistA <- dhist(locations = c(1, 3, 4), masses = c(2, 1, 1))
  dhistB <- dhist(locations = c(0, 2, 4, 5), masses = c(0.5, 2, 0.5, 1))
  
  # Set up smoothed and unsmoothed versions of histograms
  smoothing_window_width <- 1
  dhistA_unsmoothed <- as_unsmoothed_dhist(dhistA)
  dhistB_unsmoothed <- as_unsmoothed_dhist(dhistB)
  dhistA_smoothed <- as_smoothed_dhist(dhistA, smoothing_window_width)
  dhistB_smoothed <- as_smoothed_dhist(dhistB, smoothing_window_width)
  
  # Set expected area
  expected_area_unsmoothed <- 4
  expected_area_smoothed <- 3
  
  # Generate ecmfs
  ecmfA_unsmoothed <- dhist_ecmf(dhistA_unsmoothed)
  ecmfB_unsmoothed <- dhist_ecmf(dhistB_unsmoothed)
  ecmfA_smoothed <- dhist_ecmf(dhistA_smoothed)
  ecmfB_smoothed <- dhist_ecmf(dhistB_smoothed)

  # Calculate area between ECMFs
  actual_area_unsmoothed <- area_between_dhist_ecmfs(ecmfA_unsmoothed, ecmfB_unsmoothed)
  actual_area_smoothed <- area_between_dhist_ecmfs(ecmfA_smoothed, ecmfB_smoothed)
  
  # Compare caculated areas with expected areas
  expect_equal(actual_area_unsmoothed, expected_area_unsmoothed)
  expect_equal(actual_area_smoothed, expected_area_smoothed)
})

context("dhist: Area between ECMFs (non-integer normalised dhists)")
test_that("area_between_dhist_ecmfs returns correct value for non-integer normalised dhists", {
  
  # Previous simple integer grid where both histograms have been separately 
  # normalised to unit mass and variance. Has locations and masses at a range 
  # of floating point locations. Has bowties, triangles and trapeziums.
  dhistA <- dhist(locations = c(1, 3, 4), masses = c(2, 1, 1))
  dhistB <- dhist(locations = c(0, 2, 4, 5), masses = c(0.5, 2, 0.5, 1))

  dhistA <- normalise_dhist_mass(normalise_dhist_variance(dhistA))
  dhistB <- normalise_dhist_mass(normalise_dhist_variance(dhistB))
  
  # Set up smoothed and unsmoothed versions of histograms
  smoothing_window_width <- 1
  dhistA_unsmoothed <- as_unsmoothed_dhist(dhistA)
  dhistB_unsmoothed <- as_unsmoothed_dhist(dhistB)
  dhistA_smoothed <- as_smoothed_dhist(dhistA, smoothing_window_width)
  dhistB_smoothed <- as_smoothed_dhist(dhistB, smoothing_window_width)
  
  # Generate ecmfs
  ecmfA_unsmoothed <- dhist_ecmf(dhistA_unsmoothed)
  ecmfB_unsmoothed <- dhist_ecmf(dhistB_unsmoothed)
  ecmfA_smoothed <- dhist_ecmf(dhistA_smoothed)
  ecmfB_smoothed <- dhist_ecmf(dhistB_smoothed)
  
  # Define some functions to make calculation of manually measured areas easier
  rectangle_area <- function(width, height) {
    return(width * height)
  }
  triangle_area <- function(base, height) {
    return(0.5 * base * height)
  }
  trapezium_area <- function(side_a, side_b, height) {
    return(0.5 * (side_a + side_b) * height)
  }
  # Measurements of expected area between ECMFs done by hand by printing
  # normalised ECMFs on a grid with x-spacing of 0.02 and y-spacing of 0.01)
  # Actual grid counts preserved in data to facilitate less tedious manual
  # checking if required
  # --- Unsmoothed ---
  area_A_unsmoothed <- rectangle_area(width = 10*0.02, height = 12.5*0.01)
  area_B_unsmoothed <- rectangle_area(width = 50.5*0.02, height = 37.5*0.01)
  area_C_unsmoothed <- rectangle_area(width = 26*0.02, height = 12.5*0.01)
  area_D_unsmoothed <- rectangle_area(width = 34.5*0.02, height = 12.5*0.01)
  area_E_unsmoothed <- rectangle_area(width = 26.5*0.02, height = 25*0.01)
  expected_area_unsmoothed <- 
    sum(area_A_unsmoothed, area_B_unsmoothed, area_C_unsmoothed, 
        area_D_unsmoothed, area_E_unsmoothed)
  # --- Smoothed ---
  area_A_smoothed <- triangle_area(base = 2.75*0.01, height = 6.5*0.02)
  area_B_smoothed <- triangle_area(base = 2.75*0.01, height = 3*0.02)
  area_C_smoothed <- triangle_area(base = 18.5*0.01, height = 21*0.02)
  area_D_smoothed <- trapezium_area(side_a = 18.5*0.01, side_b = 37.5*0.01, height = 14.5*0.02)
  area_E_smoothed <- trapezium_area(side_a = 37.5*0.01, side_b = 37.5*0.01, height = 16*0.02)
  area_F_smoothed <- triangle_area(base = 37.5*0.01, height = 22.5*0.02)
  area_G_smoothed <- triangle_area(base = 7.5*0.01, height = 8*0.02)
  area_H_smoothed <- triangle_area(base = 7.5*0.01, height = 11*0.02)
  area_I_smoothed <- triangle_area(base = 12.5*0.01, height = 19.5*0.02)
  area_J_smoothed <- trapezium_area(side_a = 12.5*0.01, side_b = 20*0.01, height = 30.5*0.02)
  area_K_smoothed <- trapezium_area(side_a = 20*0.01, side_b = 18*0.01, height = 8*0.02)
  area_L_smoothed <- triangle_area(base = 18*0.01, height = 22*0.02)
  expected_area_smoothed <- 
    sum(area_A_smoothed, area_B_smoothed, area_C_smoothed, area_D_smoothed, 
        area_E_smoothed, area_F_smoothed, area_G_smoothed, area_H_smoothed,
        area_I_smoothed, area_J_smoothed, area_K_smoothed, area_L_smoothed)
  
  # Calculate area between ECMFs
  actual_area_unsmoothed <- area_between_dhist_ecmfs(ecmfA_unsmoothed, ecmfB_unsmoothed)
  actual_area_smoothed <- area_between_dhist_ecmfs(ecmfA_smoothed, ecmfB_smoothed)
  
  # Compare caculated areas with expected areas
  expect_equalish_manual <- function(actual, expected, relative_tolerance) {
    relative_diff <- abs(actual - expected) / expected
    expect_lte(relative_diff, relative_tolerance)
  }
  
  # Given manual measurement of areas between curves, consider area correct
  # if actual and expected areas are within 1% of each other
  expect_equalish_manual(actual_area_unsmoothed, expected_area_unsmoothed, 0.01)
  expect_equalish_manual(actual_area_smoothed, expected_area_smoothed, 0.01)
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
