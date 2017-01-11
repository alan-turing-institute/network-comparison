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

test_that("discrete_hist generates correct discrete histograms", {
  # Method for generating random observations containing specific locations a 
  # specific number of times
  random_observations <- function(locations, counts) {
    # Construct vector containing each location replicated "count" times
    observations <- purrr::simplify(purrr::map2(locations, counts, rep))
    # Randomise the order of the observations
    sample(observations)
  }
  
  # Set parameters for generation of random observation sets
  num_observations <- 15
  location_range <- -50:50
  count_range <- 0:100
  
  # Generate random observation sets
  locations <- sample(location_range, num_observations, replace = FALSE)
  counts <- sample(count_range, num_observations, replace = TRUE)
  observations <- random_observations(locations, counts)
  
  # Generate discrete histograms
  hist <- discrete_hist(observations)

  # Check that histogram locations and counts match those used to generate the 
  # observations
  sort_out <- sort(locations, index.return = TRUE)
  sorted_locations <- sort_out$x
  sorted_counts <- counts[sort_out$ix]
  
  expect_true(all.equal(hist$locations, sorted_locations))
  expect_true(all.equal(hist$counts, sorted_counts))
  
})