library("purrr")
context("dhist:Discrete histogram")

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
