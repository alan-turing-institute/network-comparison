
context("Measures Netdis: Test adaptive binning")
test_that("adaptive_breaks merges 2 lowest bins where only first bin is below minimum", {
  min_count <- 5
  x <- c(1.5, rep(2.2, min_count), rep(3.5, min_count), rep(4.5, min_count), 
         rep(5.5, min_count), rep(6.5, min_count + 1))
  initial_breaks <- 1:7
  final_breaks_actual <- adaptive_breaks(x, min_count = min_count, breaks = initial_breaks)
  final_breaks_expected <- c(1, 3, 4, 5, 6, 7)

  expect_equal(final_breaks_actual, final_breaks_expected)
})

test_that("adaptive_breaks merges 3 lowest bins where lowest 2 combined are below minimum", {
  min_count <- 5
  x <- c(1.5, rep(2.2, 2), rep(3.5, min_count), rep(4.5, min_count), 
         rep(5.5, min_count), rep(6.5, min_count + 1))
  initial_breaks <- 1:7
  final_breaks_actual <- adaptive_breaks(x, min_count = min_count, breaks = initial_breaks)
  final_breaks_expected <- c(1, 4, 5, 6, 7)
  
  expect_equal(final_breaks_actual, final_breaks_expected)
})

test_that("adaptive_breaks merges pair of bins in middle", {
  min_count <- 5
  x <- c(rep(1.6, min_count), rep(2.2, min_count), rep(3.5, 2), rep(4.5, 3), 
         rep(5.5, min_count), rep(6.5, min_count + 1))
  initial_breaks <- 1:7
  final_breaks_actual <- adaptive_breaks(x, min_count = min_count, breaks = initial_breaks)
  final_breaks_expected <- c(1, 2, 3, 5, 6, 7)
  
  expect_equal(final_breaks_actual, final_breaks_expected)
})

test_that("adaptive_breaks merges two spearated pairs of bins in middle", {
  min_count <- 5
  x <- c(rep(1.6, min_count), rep(2.2, 2), rep(3.5, 3), rep(4.5, min_count), 
         rep(5.5, 3), rep(6.5, 2), rep(7.8, min_count))
  initial_breaks <- 1:8
  final_breaks_actual <- adaptive_breaks(x, min_count = min_count, breaks = initial_breaks)
  final_breaks_expected <- c(1, 2, 4, 5, 7, 8)
  
  expect_equal(final_breaks_actual, final_breaks_expected)
})

test_that("adaptive_breaks merges 2 uppermost bins where both are below minimum", {
  min_count <- 5
  x <- c(rep(1.5, min_count), rep(2.2, min_count), rep(3.5, min_count),
         rep(4.5, min_count), rep(5.5, 2), rep(6.5, 3))
  initial_breaks <- 1:7
  final_breaks_actual <- adaptive_breaks(x, min_count = min_count, breaks = initial_breaks)
  final_breaks_expected <- c(1, 2,3, 4, 5, 7)

  expect_equal(final_breaks_actual, final_breaks_expected)
})

test_that("adaptive_breaks merges 2 uppermost bins where only last bin is below minimum", {
  min_count <- 5
  x <- c(rep(1.5, min_count), rep(2.2, min_count), rep(3.5, min_count),
         rep(4.5, min_count), rep(5.5, min_count), rep(6.5, 3))
  initial_breaks <- 1:7
  final_breaks_actual <- adaptive_breaks(x, min_count = min_count, breaks = initial_breaks)
  final_breaks_expected <- c(1, 2,3, 4, 5, 7)
  
  expect_equal(final_breaks_actual, final_breaks_expected)
})