
self_net_emd <- function(histogram, shift, method) {
  netemd_one_to_one(dhists_1 = histogram, dhists_2 = shift_dhist(histogram, shift), method = method)
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

expect_self_netemd_correct <- function(histogram, shift, method,
                                        return_details = FALSE) {
  self_net_emd <- netemd_one_to_one(dhists_1 = histogram, dhists_2 = shift_dhist(histogram, shift),
    method = method, return_details = return_details
  )
  loc <- histogram$locations
  mass <- histogram$masses
  var <- sum(loc * loc * mass) / sum(mass) - (sum(loc * mass) / sum(mass))^2
  expected <- list(
    net_emd = 0, min_emds = 0, min_offsets = shift,
    min_offsets_std = 0
  )
  expect_equal(self_net_emd, expected)
}

locations <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)
masses <- c(0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0)
histogram <- dhist(locations = locations, masses = masses)
expect_self_netemd_correct(histogram,
  shift = 1, "optimise",
  return_details = TRUE
)
expect_self_netemd_correct(histogram,
  shift = 1, "exhaustive",
  return_details = TRUE
)
expect_self_netemd_correct(histogram,
  shift = 0.5, "optimise",
  return_details = TRUE
)
expect_self_netemd_correct(histogram,
  shift = 0.5, "exhaustive",
  return_details = TRUE
)
expect_self_netemd_correct(histogram,
  shift = 0.1, "optimise",
  return_details = TRUE
)
expect_self_netemd_correct(histogram,
  shift = 0.1, "exhaustive",
  return_details = TRUE
)
expect_self_netemd_correct(histogram,
  shift = 0.05, "optimise",
  return_details = TRUE
)
expect_self_netemd_correct(histogram,
  shift = 0.05, "exhaustive",
  return_details = TRUE
)
expect_self_netemd_correct(histogram,
  shift = 0.01, "optimise",
  return_details = TRUE
)
expect_self_netemd_correct(histogram,
  shift = 0.01, "exhaustive",
  return_details = TRUE
)
expect_self_netemd_correct(histogram,
  shift = 0, "optimise",
  return_details = TRUE
)
expect_self_netemd_correct(histogram,
  shift = 0, "exhaustive",
  return_details = TRUE
)

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
    netemd_one_to_one(dhists_1 = dhist, dhists_2 = dhist, method = "optimise")
  })
  purrr::walk(actuals_opt, function(actual) {
    expect_equal(actual, expected)
  })

  actuals_exhaustive_default <- purrr::map(rand_dhists, function(dhist) {
    netemd_one_to_one(dhists_1 = dhist, dhists_2 = dhist, method = "exhaustive")
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

  netemd_offset_self <- function(dhist, offsets, method) {
    netemds <- purrr::map_dbl(offsets, function(offset) {
      netemd_one_to_one(dhists_1 = dhist, dhists_2 = shift_dhist(dhist, offset), method = method)
    })
    return(netemds)
  }

  expected <- 0
  actuals_list_opt <- purrr::map2(
    rand_dhists, offset_lists,
    function(dhist, offsets) {
      netemd_offset_self(dhist, offsets, method = "optimise")
    }
  )
  purrr::walk(actuals_list_opt, function(actuals) {
    purrr::walk(actuals, function(actual) {
      expect_equal(actual, expected)
    })
  })
  actuals_list_exhaustive <- purrr::map2(
    rand_dhists, offset_lists,
    function(dhist, offsets) {
      netemd_offset_self(dhist, offsets, method = "exhaustive")
    }
  )
  purrr::walk(actuals_list_exhaustive, function(actuals) {
    purrr::walk(actuals, function(actual) {
      expect_equal(actual, expected)
    })
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
    return(seq(mu - 5 * sigma, mu + 5 * sigma, length.out = num_bins))
  }

  rand_dhists <- purrr::map2(mus, sigmas, function(mu, sigma) {
    locations <- rand_locations(mu, sigma)
    masses <- dnorm(locations, mean = mu, sd = sigma)
    return(dhist(masses = masses, locations = locations))
  })

  offset_lists <- replicate(num_hists, offsets, simplify = FALSE)

  expect_self_netemd_correct <-
    function(histogram, shift, method, return_details = FALSE) {
      self_net_emd <- netemd_one_to_one(dhists_1 = histogram, dhists_2 = shift_dhist(histogram, shift),method = method, return_details = return_details
      )
      loc <- histogram$locations
      mass <- histogram$masses
      var <- sum(loc * loc * mass) / sum(mass) - (sum(loc * mass) / sum(mass))^2
      expected <- list(
        net_emd = 0, min_emds = 0, min_offsets = shift,
        min_offsets_std = 0
      )
      expect_equal(self_net_emd, expected)
    }

  purrr::walk2(rand_dhists, offset_lists, function(dhist, offsets) {
    purrr::walk(offsets, function(offset) {
      expect_self_netemd_correct(dhist, offset, "optimise",
        return_details = TRUE
      )
    })
  })

  purrr::walk2(rand_dhists, offset_lists, function(dhist, offsets) {
    purrr::walk(offsets, function(offset) {
      expect_self_netemd_correct(dhist, offset, "exhaustive",
        return_details = TRUE
      )
    })
  })
})

test_that("net_emd returns analytically derived non-zero solutions for distributions
          where the analytical solution is known", {
  # Helper functions to create dhists for a given value of "p"
  two_bin_dhist <- function(p) {
    dhist(locations = c(0, 1), masses = c(p, 1 - p))
  }
  three_bin_dhist <- function(p) {
    dhist(locations = c(-1, 0, 1), masses = c(0.5 * p * (1 - p), 1 - (p * (1 - p)), 0.5 * p * (1 - p)))
  }

  # Helper function to test actual vs expected
  test_pair <- function(p, expected) {
    dhistA <- two_bin_dhist(p)
    dhistB <- three_bin_dhist(p)
    expect_equal(netemd_one_to_one(dhists_1 = dhistA, dhists_2 = dhistB, method = "exhaustive"), expected, tolerance = 1e-12)
    # Even setting the stats::optimise method tolerance to machine double precision, the
    # optimised NetEMD is ~1e-09, so set a slightly looser tolerance here
    expect_equal(netemd_one_to_one(dhists_1 = dhistA, dhists_2 = dhistB, method = "optimise"), expected, tolerance = 1e-08)
  }

  # Test for p values with analytically calculated NetEMD
  test_pair(1 / 2, 1)
  test_pair(1 / 3, 1 / sqrt(2))
  test_pair(1 / 5, 1 / 2)
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
      expect_equalish(netemd_one_to_one(dhists_1 = gdd_Ox, dhists_2 = gdd_Ox,
        method = "optimise",
        smoothing_window_width = 0
      ), 0)
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
    format = "ncol", pattern = "*"
  )
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
    format = "ncol", pattern = "*"
  )
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
      expect_equalish(netemd_one_to_one(dhists_1 = gdd_Ox, dhists_2 = gdd_Ox,
        method = "optimise",
        smoothing_window_width = 0
      ), 0)
    })
  })
})

context("Measures NetEMD: All graphs in directory")
test_that("netemd_many_to_many works", {
  # Set source directory and file properties for Virus PPI graph edge files
  source_dir <- system.file(file.path("extdata", "VRPINS"), package = "netdist")
  edge_format <- "ncol"
  file_pattern <- ".txt"

  # Set number of threads to use at once for parallel processing.
  num_threads <- getOption("mc.cores", 2L)

  # Use previously tested GDD code to generate inputs to function under test
  gdds_orbits_g4 <- gdd_for_all_graphs(
    source_dir = source_dir, format = edge_format, pattern = file_pattern,
    feature_type = "orbit", max_graphlet_size = 4
  )
  gdds_orbits_g5 <- gdd_for_all_graphs(
    source_dir = source_dir, format = edge_format, pattern = file_pattern,
    feature_type = "orbit", max_graphlet_size = 5
  )
  gdds_graphlets_g4 <- gdd_for_all_graphs(
    source_dir = source_dir, format = edge_format, pattern = file_pattern,
    feature_type = "graphlet", max_graphlet_size = 4
  )
  gdds_graphlets_g5 <- gdd_for_all_graphs(
    source_dir = source_dir, format = edge_format, pattern = file_pattern,
    feature_type = "graphlet", max_graphlet_size = 5
  )
  gdds_graphlets_g4_e1 <- gdd_for_all_graphs(
    source_dir = source_dir, format = edge_format, pattern = file_pattern,
    feature_type = "graphlet", max_graphlet_size = 4, ego_neighbourhood_size = 1
  )
  gdds_graphlets_g5_e1 <- gdd_for_all_graphs(
    source_dir = source_dir, format = edge_format, pattern = file_pattern,
    feature_type = "graphlet", max_graphlet_size = 5, ego_neighbourhood_size = 1
  )
  gdds_graphlets_g4_e2 <- gdd_for_all_graphs(
    source_dir = source_dir, format = edge_format, pattern = file_pattern,
    feature_type = "graphlet", max_graphlet_size = 4, ego_neighbourhood_size = 2
  )
  gdds_graphlets_g5_e2 <- gdd_for_all_graphs(
    source_dir = source_dir, format = edge_format, pattern = file_pattern,
    feature_type = "graphlet", max_graphlet_size = 5, ego_neighbourhood_size = 2
  )

  # Use previously tested NetEMD function to generate expected NetEMD scores
  # individually and combine into expected output for code under test
  expected_netemd_fn <- function(gdds) {
    list(
      netemds = c(
        netemd_one_to_one(dhists_1 = gdds$EBV, dhists_2 =  gdds$ECL), netemd_one_to_one(dhists_1 =gdds$EBV, dhists_2 = gdds$HSV),
        netemd_one_to_one(dhists_1 = gdds$EBV, dhists_2 = gdds$KSHV), netemd_one_to_one(dhists_1 =gdds$EBV, dhists_2 = gdds$VZV),
        netemd_one_to_one(dhists_1 = gdds$ECL, dhists_2 = gdds$HSV), netemd_one_to_one(dhists_1 =gdds$ECL, dhists_2 = gdds$KSHV),
        netemd_one_to_one(dhists_1 = gdds$ECL, dhists_2 = gdds$VZV), netemd_one_to_one(dhists_1 =gdds$HSV, dhists_2 = gdds$KSHV),
        netemd_one_to_one(dhists_1 = gdds$HSV, dhists_2 = gdds$VZV), netemd_one_to_one(dhists_1 =gdds$KSHV, dhists_2 = gdds$VZV)
      ),
      comp_spec = cross_comparison_spec(gdds)
    )
  }

  # Comparison function for clarity
  compare_fn <- function(gdds) {
    expect_equal(netemd_many_to_many(dhists=gdds), expected_netemd_fn(gdds))
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
