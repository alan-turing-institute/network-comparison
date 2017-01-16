library("purrr")
context("Utility Functions")

test_that("rotl_vec rotates left by specified number of places", {
  test_vec <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
  expect_equal(rotl_vec(test_vec, -13), c(80, 90, 100, 10, 20, 30, 40, 50, 60, 70))
  expect_equal(rotl_vec(test_vec, -12), c(90, 100, 10, 20, 30, 40, 50, 60, 70, 80))
  expect_equal(rotl_vec(test_vec, -11), c(100, 10, 20, 30, 40, 50, 60, 70, 80, 90))
  expect_equal(rotl_vec(test_vec, -10), c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100))
  expect_equal(rotl_vec(test_vec, -9), c(20, 30, 40, 50, 60, 70, 80, 90, 100, 10))
  expect_equal(rotl_vec(test_vec, -8), c(30, 40, 50, 60, 70, 80, 90, 100, 10, 20))
  expect_equal(rotl_vec(test_vec, -7), c(40, 50, 60, 70, 80, 90, 100, 10, 20, 30))
  expect_equal(rotl_vec(test_vec, -6), c(50, 60, 70, 80, 90, 100, 10, 20, 30, 40))
  expect_equal(rotl_vec(test_vec, -5), c(60, 70, 80, 90, 100, 10, 20, 30, 40, 50))
  expect_equal(rotl_vec(test_vec, -4), c(70, 80, 90, 100, 10, 20, 30, 40, 50, 60))
  expect_equal(rotl_vec(test_vec, -3), c(80, 90, 100, 10, 20, 30, 40, 50, 60, 70))
  expect_equal(rotl_vec(test_vec, -2), c(90, 100, 10, 20, 30, 40, 50, 60, 70, 80))
  expect_equal(rotl_vec(test_vec, -1), c(100, 10, 20, 30, 40, 50, 60, 70, 80, 90))
  expect_equal(rotl_vec(test_vec, 0), c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100))
  expect_equal(rotl_vec(test_vec, 1), c(20, 30, 40, 50, 60, 70, 80, 90, 100, 10))
  expect_equal(rotl_vec(test_vec, 2), c(30, 40, 50, 60, 70, 80, 90, 100, 10, 20))
  expect_equal(rotl_vec(test_vec, 3), c(40, 50, 60, 70, 80, 90, 100, 10, 20, 30))
  expect_equal(rotl_vec(test_vec, 4), c(50, 60, 70, 80, 90, 100, 10, 20, 30, 40))
  expect_equal(rotl_vec(test_vec, 5), c(60, 70, 80, 90, 100, 10, 20, 30, 40, 50))
  expect_equal(rotl_vec(test_vec, 6), c(70, 80, 90, 100, 10, 20, 30, 40, 50, 60))
  expect_equal(rotl_vec(test_vec, 7), c(80, 90, 100, 10, 20, 30, 40, 50, 60, 70))
  expect_equal(rotl_vec(test_vec, 8), c(90, 100, 10, 20, 30, 40, 50, 60, 70, 80))
  expect_equal(rotl_vec(test_vec, 9), c(100, 10, 20, 30, 40, 50, 60, 70, 80, 90))
  expect_equal(rotl_vec(test_vec, 10), c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100))
  expect_equal(rotl_vec(test_vec, 11), c(20, 30, 40, 50, 60, 70, 80, 90, 100, 10))
  expect_equal(rotl_vec(test_vec, 12), c(30, 40, 50, 60, 70, 80, 90, 100, 10, 20))
  expect_equal(rotl_vec(test_vec, 13), c(40, 50, 60, 70, 80, 90, 100, 10, 20, 30))
})


test_that("rotr_vec rotates right by specified number of places", {
  test_vec <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
  expect_equal(rotr_vec(test_vec, 13), c(80, 90, 100, 10, 20, 30, 40, 50, 60, 70))
  expect_equal(rotr_vec(test_vec, 12), c(90, 100, 10, 20, 30, 40, 50, 60, 70, 80))
  expect_equal(rotr_vec(test_vec, 11), c(100, 10, 20, 30, 40, 50, 60, 70, 80, 90))
  expect_equal(rotr_vec(test_vec, 10), c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100))
  expect_equal(rotr_vec(test_vec, 9), c(20, 30, 40, 50, 60, 70, 80, 90, 100, 10))
  expect_equal(rotr_vec(test_vec, 8), c(30, 40, 50, 60, 70, 80, 90, 100, 10, 20))
  expect_equal(rotr_vec(test_vec, 7), c(40, 50, 60, 70, 80, 90, 100, 10, 20, 30))
  expect_equal(rotr_vec(test_vec, 6), c(50, 60, 70, 80, 90, 100, 10, 20, 30, 40))
  expect_equal(rotr_vec(test_vec, 5), c(60, 70, 80, 90, 100, 10, 20, 30, 40, 50))
  expect_equal(rotr_vec(test_vec, 4), c(70, 80, 90, 100, 10, 20, 30, 40, 50, 60))
  expect_equal(rotr_vec(test_vec, 3), c(80, 90, 100, 10, 20, 30, 40, 50, 60, 70))
  expect_equal(rotr_vec(test_vec, 2), c(90, 100, 10, 20, 30, 40, 50, 60, 70, 80))
  expect_equal(rotr_vec(test_vec, 1), c(100, 10, 20, 30, 40, 50, 60, 70, 80, 90))
  expect_equal(rotr_vec(test_vec, 0), c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100))
  expect_equal(rotr_vec(test_vec, -1), c(20, 30, 40, 50, 60, 70, 80, 90, 100, 10))
  expect_equal(rotr_vec(test_vec, -2), c(30, 40, 50, 60, 70, 80, 90, 100, 10, 20))
  expect_equal(rotr_vec(test_vec, -3), c(40, 50, 60, 70, 80, 90, 100, 10, 20, 30))
  expect_equal(rotr_vec(test_vec, -4), c(50, 60, 70, 80, 90, 100, 10, 20, 30, 40))
  expect_equal(rotr_vec(test_vec, -5), c(60, 70, 80, 90, 100, 10, 20, 30, 40, 50))
  expect_equal(rotr_vec(test_vec, -6), c(70, 80, 90, 100, 10, 20, 30, 40, 50, 60))
  expect_equal(rotr_vec(test_vec, -7), c(80, 90, 100, 10, 20, 30, 40, 50, 60, 70))
  expect_equal(rotr_vec(test_vec, -8), c(90, 100, 10, 20, 30, 40, 50, 60, 70, 80))
  expect_equal(rotr_vec(test_vec, -9), c(100, 10, 20, 30, 40, 50, 60, 70, 80, 90))
  expect_equal(rotr_vec(test_vec, -10), c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100))
  expect_equal(rotr_vec(test_vec, -11), c(20, 30, 40, 50, 60, 70, 80, 90, 100, 10))
  expect_equal(rotr_vec(test_vec, -12), c(30, 40, 50, 60, 70, 80, 90, 100, 10, 20))
  expect_equal(rotr_vec(test_vec, -13), c(40, 50, 60, 70, 80, 90, 100, 10, 20, 30))
})

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
    hist <- discrete_hist(observations)
    
    # discrete_hist will drop bins with zero counts, so remove these from the 
    # expected data (not necessary now we've restricted counts to be >= 1, but 
    # the bug where we generated test locations with zero counts was so annoying
    # to identify that we're going with a belt and braces approach)
    non_zero_count_indexes <- counts != 0
    expected_locations <- locations[non_zero_count_indexes]
    expected_counts <- counts[non_zero_count_indexes]
    # discrete_hist will return results with bins ordered by ascending location, 
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
