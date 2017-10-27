context("Measures NetEMD: NetEMD")
# NetEMD: Property-based tests
test_that("net_emd returns 0 when comparing an integer location histogram offset
          against itself", {

  self_net_emd <- function(histogram, shift, method) {
    net_emd(histogram, shift_dhist(histogram, shift), method = method)
  }
  expected <- 0

  locations <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)
  masses <- c(0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0)
  histogram <- dhist(locations = locations, masses = masses)

  expect_equal(self_net_emd(histogram, shift = 1, "optimise"), expected)
  expect_equal(self_net_emd(histogram, shift = 1, "exhaustive"), expected)
  expect_equal(self_net_emd(histogram, shift = 0.5, "optimise"), expected)
  expect_equal(self_net_emd(histogram, shift = 0.5, "exhaustive"), expected)
  expect_equal(self_net_emd(histogram, shift = 0.1, "optimise"), expected)
  expect_equal(self_net_emd(histogram, shift = 0.1, "exhaustive"), expected)
  expect_equal(self_net_emd(histogram, shift = 0.05, "optimise"), expected)
  expect_equal(self_net_emd(histogram, shift = 0.05, "exhaustive"), expected)
  expect_equal(self_net_emd(histogram, shift = 0.01, "optimise"), expected)
  expect_equal(self_net_emd(histogram, shift = 0.01, "exhaustive"), expected)
  expect_equal(self_net_emd(histogram, shift = 0, "optimise"), expected)
  expect_equal(self_net_emd(histogram, shift = 0, "exhaustive"), expected)
})

test_that("net_emd returns min_emd = 0 and min_offset = 0 when comparing an
          integer location histogram offset against itself", {

  expect_self_net_emd_correct <- function(histogram, shift, method,
                                          return_details = FALSE) {
    self_net_emd <- net_emd(histogram, shift_dhist(histogram, shift),
                            method = method, return_details = return_details)
    expected <- list(net_emd = 0, min_emds = 0, min_offsets = shift)
    expect_equal(self_net_emd, expected)
  }

  locations <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)
  masses <- c(0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0)
  histogram <- dhist(locations = locations, masses = masses)

  expect_self_net_emd_correct(histogram, shift = 1, "optimise",
                              return_details = TRUE)
  expect_self_net_emd_correct(histogram, shift = 1, "exhaustive",
                              return_details = TRUE)
  expect_self_net_emd_correct(histogram, shift = 0.5, "optimise",
                              return_details = TRUE)
  expect_self_net_emd_correct(histogram, shift = 0.5, "exhaustive",
                              return_details = TRUE)
  expect_self_net_emd_correct(histogram, shift = 0.1, "optimise",
                              return_details = TRUE)
  expect_self_net_emd_correct(histogram, shift = 0.1, "exhaustive",
                              return_details = TRUE)
  expect_self_net_emd_correct(histogram, shift = 0.05, "optimise",
                              return_details = TRUE)
  expect_self_net_emd_correct(histogram, shift = 0.05, "exhaustive",
                              return_details = TRUE)
  expect_self_net_emd_correct(histogram, shift = 0.01, "optimise",
                              return_details = TRUE)
  expect_self_net_emd_correct(histogram, shift = 0.01, "exhaustive",
                              return_details = TRUE)
  expect_self_net_emd_correct(histogram, shift = 0, "optimise",
                              return_details = TRUE)
  expect_self_net_emd_correct(histogram, shift = 0, "exhaustive",
                              return_details = TRUE)
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
  purrr::walk(actuals_opt, function(actual) {
    expect_equal(actual, expected)
    })

  actuals_exhaustive_default <- purrr::map(rand_dhists, function(dhist) {
    net_emd(dhist, dhist, method = "exhaustive")
    })
  purrr::walk(actuals_exhaustive_default, function(actual) {
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
  purrr::walk(actuals_list_opt, function(actuals) {
        purrr::walk(actuals, function(actual) {
          expect_equal(actual, expected)})
  })
  actuals_list_exhaustive <- purrr::map2(rand_dhists, offset_lists,
                                   function(dhist, offsets) {
    net_emd_offset_self(dhist, offsets, method = "exhaustive")})
  purrr::walk(actuals_list_exhaustive, function(actuals) {
    purrr::walk(actuals, function(actual) {expect_equal(actual, expected)})
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
    function(histogram, shift, method, return_details = FALSE) {
    self_net_emd <- net_emd(histogram, shift_dhist(histogram, shift),
                            method, return_details)
    expected <- list(net_emd = 0, min_emds = 0, min_offsets = shift)
    expect_equal(self_net_emd, expected)
  }

  purrr::walk2(rand_dhists, offset_lists, function(dhist, offsets) {
    purrr::walk(offsets, function(offset){
      expect_self_net_emd_correct(dhist, offset, "optimise",
                                  return_details = TRUE)
    })
  })

  purrr::walk2(rand_dhists, offset_lists, function(dhist, offsets) {
    purrr::walk(offsets, function(offset){
      expect_self_net_emd_correct(dhist, offset, "exhaustive",
                                  return_details = TRUE)
    })
  })
})

test_that("net_emd returns analytically derived non-zero solutions for distributions
          where the analytical solution is known", {
  # Helper functions to create dhists for a given value of "p"
  two_bin_dhist <- function(p) {
    dhist(locations = c(0, 1), masses = c(p, 1-p))
  }
  three_bin_dhist <- function(p) {
  dhist(locations = c(-1, 0, 1), masses = c(0.5*p*(1-p), 1-(p*(1-p)), 0.5*p*(1-p)))
  }
  
  # Helper function to test actual vs expected
  test_pair <- function(p, expected) {
    dhistA <- two_bin_dhist(p)
    dhistB <- three_bin_dhist(p)
    expect_equal(net_emd(dhistA, dhistB, method = "exhaustive"), expected, tolerance = 1e-12)
    # Even setting the stats::optimise method tolerance to machine double precision, the
    # optimised NetEMD is ~1e-09, so set a slightly looser tolerance here
    expect_equal(net_emd(dhistA, dhistB, method = "optimise"), expected, tolerance = 1e-08)
  }
  
  # Test for p values with analytically calculated NetEMD
  test_pair(1/2, 1)
  test_pair(1/3, 1/sqrt(2))
  test_pair(1/5, 1/2)
   
})

test_that("net_emd gives same results as underlying min_emd_* methods", {
  
  bin_masses1 <- c(2.124749, 3.453120, 5.936767, 6.228690, 8.057908, 10.463148, 11.394959, 11.602371, 11.922119)
  bin_centres1 <- c(-3.876220, -3.731132, -3.630938, -3.096881, -1.436239, -0.729344, -0.226891, 3.570451, 3.874552)
  bin_masses2 <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  bin_centres2 <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)
  
  # NetEMD normalises histograms to unit mas and variance, so we must do this in advance if we are to expect
  # the min_emd_* methods to give the same results
  dhist1 <- normalise_dhist_variance(normalise_dhist_mass(dhist(masses = bin_masses1, locations = bin_centres1)))
  dhist2 <- normalise_dhist_variance(normalise_dhist_mass(dhist(masses = bin_masses2, locations = bin_centres2)))
  
  # Helper function for comparing NetEMD and MinEMD output
  expect_equivalent <- function(net_emd, min_emd) {
    expect_equal(length(net_emd$min_emds), 1)
    expect_equal(length(net_emd$min_offsets), 1)
    expect_equal(net_emd$min_emds[[1]], min_emd$min_emd)
    expect_equal(net_emd$min_offsets[[1]], min_emd$min_offset)
  }
  
  # No method specific argments (use defaults from specific methods)
  expect_equivalent(net_emd(dhist1, dhist2, method = "exhaustive", return_details = TRUE), min_emd_exhaustive(dhist1, dhist2))
  expect_equivalent(net_emd(dhist1, dhist2, method = "optimise", return_details = TRUE), min_emd_optimise(dhist1, dhist2))
  expect_equivalent(net_emd(dhist1, dhist2, method = "fixed_step", return_details = TRUE), min_emd_fixed_step(dhist1, dhist2))
  
  # Providing all method-specific arguments
  # Exhaustive: No method-specific arguments to provide
  # Optimise: Vary "tolerance" parameter
  expect_equivalent(net_emd(dhist1, dhist2, method = "optimise", tolerance = 0.1, return_details = TRUE), min_emd_optimise(dhist1, dhist2, tolerance = 0.1))
  expect_equivalent(net_emd(dhist1, dhist2, method = "optimise", tolerance = 0.01, return_details = TRUE), min_emd_optimise(dhist1, dhist2, tolerance = 0.01))
  expect_equivalent(net_emd(dhist1, dhist2, method = "optimise", tolerance = 0.001, return_details = TRUE), min_emd_optimise(dhist1, dhist2, tolerance = 0.001))
  # Fixed step: Vary "step size" parameter
  expect_equivalent(net_emd(dhist1, dhist2, method = "fixed_step", step_size = 0.1, return_details = TRUE), min_emd_fixed_step(dhist1, dhist2, step_size = 0.1))
  expect_equivalent(net_emd(dhist1, dhist2, method = "fixed_step", step_size = 0.01, return_details = TRUE), min_emd_fixed_step(dhist1, dhist2, step_size = 0.01))
  expect_equivalent(net_emd(dhist1, dhist2, method = "fixed_step", step_size = 0.001, return_details = TRUE), min_emd_fixed_step(dhist1, dhist2, step_size = 0.001))
  
})

context("Measures NetEMD: Virus PPI (EMD)")
# EMD and NET_EMD: Virus PPI datasets
test_that("emd return 0 when comparing graphlet orbit degree distributions of
          virus PPI graphs to themselves", {
            # Load viurs PPI network data in ORCA-compatible edge list format
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
            random_graphs <- read_simple_graphs(
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
            random_graphs <- read_simple_graphs(
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


