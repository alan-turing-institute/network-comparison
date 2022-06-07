context("Graph binning: Adaptive binning")
test_that(
  "adaptive_breaks merges 2 lowest bins where only first bin is below minimum",
  {
    min_count <- 5
    x <- c(
      1.5, rep(2.2, min_count), rep(3.5, min_count), rep(4.5, min_count),
      rep(5.5, min_count), rep(6.5, min_count + 1)
    )
    initial_breaks <- 1:7
    final_breaks_actual <- adaptive_breaks(
      x,
      min_count = min_count,
      breaks = initial_breaks
    )
    final_breaks_expected <- c(1, 3, 4, 5, 6, 7)

    expect_equal(final_breaks_actual, final_breaks_expected)
  }
)

test_that(
  paste(
    "adaptive_breaks merges 3 lowest bins where lowest 2 combined are below",
    "minimum"
  ),
  {
    min_count <- 5
    x <- c(
      1.5, rep(2.2, 2), rep(3.5, min_count), rep(4.5, min_count),
      rep(5.5, min_count), rep(6.5, min_count + 1)
    )
    initial_breaks <- 1:7
    final_breaks_actual <- adaptive_breaks(
      x,
      min_count = min_count,
      breaks = initial_breaks
    )
    final_breaks_expected <- c(1, 4, 5, 6, 7)

    expect_equal(final_breaks_actual, final_breaks_expected)
  }
)

test_that("adaptive_breaks merges pair of bins in middle", {
  min_count <- 5
  x <- c(
    rep(1.6, min_count), rep(2.2, min_count), rep(3.5, 2), rep(4.5, 3),
    rep(5.5, min_count), rep(6.5, min_count + 1)
  )
  initial_breaks <- 1:7
  final_breaks_actual <- adaptive_breaks(
    x,
    min_count = min_count,
    breaks = initial_breaks
  )
  final_breaks_expected <- c(1, 2, 3, 5, 6, 7)

  expect_equal(final_breaks_actual, final_breaks_expected)
})

test_that("adaptive_breaks merges two spearated pairs of bins in middle", {
  min_count <- 5
  x <- c(
    rep(1.6, min_count), rep(2.2, 2), rep(3.5, 3), rep(4.5, min_count),
    rep(5.5, 3), rep(6.5, 2), rep(7.8, min_count)
  )
  initial_breaks <- 1:8
  final_breaks_actual <- adaptive_breaks(
    x,
    min_count = min_count,
    breaks = initial_breaks
  )
  final_breaks_expected <- c(1, 2, 4, 5, 7, 8)

  expect_equal(final_breaks_actual, final_breaks_expected)
})

test_that(
  "adaptive_breaks merges 2 uppermost bins where both are below minimum",
  {
    min_count <- 5
    x <- c(
      rep(1.5, min_count), rep(2.2, min_count), rep(3.5, min_count),
      rep(4.5, min_count), rep(5.5, 2), rep(6.5, 3)
    )
    initial_breaks <- 1:7
    final_breaks_actual <- adaptive_breaks(
      x,
      min_count = min_count,
      breaks = initial_breaks
    )
    final_breaks_expected <- c(1, 2, 3, 4, 5, 7)

    expect_equal(final_breaks_actual, final_breaks_expected)
  }
)

test_that(
  paste(
    "adaptive_breaks merges 2 uppermost bins where only last bin is below",
    "minimum"
  ),
  {
    min_count <- 5
    x <- c(
      rep(1.5, min_count), rep(2.2, min_count), rep(3.5, min_count),
      rep(4.5, min_count), rep(5.5, min_count), rep(6.5, 3)
    )
    initial_breaks <- 1:7
    final_breaks_actual <- adaptive_breaks(
      x,
      min_count = min_count,
      breaks = initial_breaks
    )
    final_breaks_expected <- c(1, 2, 3, 4, 5, 7)

    expect_equal(final_breaks_actual, final_breaks_expected)
  }
)

test_that("adaptive_breaks merges bins with no members with the next bin", {
  min_count <- 5
  x <- c(rep(1.5, min_count), rep(5.5, min_count), rep(6.5, min_count))
  initial_breaks <- 1:7
  final_breaks_actual <- adaptive_breaks(
    x,
    min_count = min_count,
    breaks = initial_breaks
  )
  final_breaks_expected <- c(1, 2, 6, 7)

  expect_equal(final_breaks_actual, final_breaks_expected)
})

test_that(
  paste(
    "adaptive_breaks merges 2 bins below minimum, plus the empty bins between",
    "them"
  ),
  {
    min_count <- 5
    x <- c(rep(1.5, min_count), rep(2.3, 1), rep(5.5, 4), rep(6.5, min_count))
    initial_breaks <- 1:7
    final_breaks_actual <- adaptive_breaks(
      x,
      min_count = min_count,
      breaks = initial_breaks
    )
    final_breaks_expected <- c(1, 2, 6, 7)

    expect_equal(final_breaks_actual, final_breaks_expected)
  }
)

context("Graph binning:  Adaptively binned densities")
test_that("binned_densities_adaptive works", {
  # Helper function
  test_binning <- function(densities, min_counts_per_interval, num_intervals,
                           breaks, expected_breaks, expected_interval_indexes) {
    # Set up expected output
    expected <- list(
      densities = densities,
      interval_indexes = expected_interval_indexes,
      breaks = expected_breaks
    )
    # Calculate actual output
    actual <- binned_densities_adaptive(
      densities,
      min_counts_per_interval = min_counts_per_interval,
      num_intervals = num_intervals
    )
    # Check actual matches expected
    expect_equal(actual, expected)
  }
  # Test 1:
  densities <- c(0, 0.099, 0.2, 0.299, 0.4, 0.49, 0.6, 0.699, 0.8, 0.899, 1.0)
  min_counts_per_interval <- 2
  num_intervals <- 100
  expected_breaks <- c(0, 0.1, 0.3, 0.5, 0.7, 1.0)
  expected_interval_indexes <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5)
  test_binning(densities,
    min_counts_per_interval = min_counts_per_interval,
    num_intervals = num_intervals, expected_breaks = expected_breaks,
    expected_interval_indexes = expected_interval_indexes
  )
  # Test 2:
  densities <- c(
    0,
    0.012,
    0.099,
    0.201,
    0.299,
    0.402,
    0.49,
    0.596,
    0.699,
    0.803,
    0.899,
    1.0
  )
  min_counts_per_interval <- 2
  num_intervals <- 100
  expected_breaks <- c(0, 0.02, 0.21, 0.41, 0.6, 0.81, 1.0)
  expected_interval_indexes <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6)
  test_binning(densities,
    min_counts_per_interval = min_counts_per_interval,
    num_intervals = num_intervals, expected_breaks = expected_breaks,
    expected_interval_indexes = expected_interval_indexes
  )
})

expected_binned_graphlet_counts <-
  function(graphs, binning_fn, max_graphlet_size) {
    binned_graphs <- binning_fn(graphs)
    ref_counts <- purrr::map(
      binned_graphs$graphs, count_graphlets_for_graph,
      max_graphlet_size
    )
    ref_counts
  }
