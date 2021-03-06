/*
 * This file uses the Catch unit testing library, alongside
 * testthat's simple bindings, to test a C++ function.
 *
 * This file should begin with `test` and be in the  `src/` folder.
 * `LinkingTo: testthat` must also be within the DESCRIPTION file.
 */
// Enable C++11
// [[Rcpp::plugins(cpp11)]]

// All test files should include the <testthat.h>
// header file.
#include "emd_fast_no_smoothing.h"
#include <math.h>
#include <vector>
#include <testthat.h>

// Helper function to test tolerance
bool within_tolerance(double actual, double expected, double tolerance) {
  if(actual > expected) {
    return ((actual - expected) <= tolerance);
  }
  else {
    return ((expected - actual) <= tolerance);
  }
}

// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.
context("emd_fast_no_smoothing") {
  
  // Test Kahan summation helper function
  // We also verify that the example data used to test the Kahan summation
  // results in an unacceptably large error when using naiive summation.
  test_that("Kahan summation works when naiive summation fails") {
    // NOTE: If this test fails when the `add_element_kahan()` function has not
    // been changed, there is a small chance that it may be that the package is
    // being compiled with the `-ffast-math` compiler flag set (or an aggressive
    // optimisation level that sets the fast math flag). If this flag is set, 
    // the Kahan compensation being optimised away.
    
    // Set up suitable test data
    // =========================
    // NOTE: It is surprisingly easy to come up with test data where the naiive
    // summation works fine, so the test to validate that the naiive sum results
    // in an unnacceptably high error is important. I think this is for the 
    // following reasons.
    // (A) If the start number is too large in comparison to the expected sum of
    //     the elements in the small number vector, then the floating point
    //     representation of the expected total will be the same as that of the
    //     start number.
    // (B) If the start number is too small in comparison to the elements in the
    //     small number vector, we don't run into the floating point 
    //     representation issue when we add each element to the running total.
    // NOTE: Appveyor tests with both 32 and 64 bit builds of R (i386 and ix64 
    //       respectively). Therefore values have been chosen that all fall
    //       within a 32-bit floating point range (including the expected total)
    //       while still causing the naiive summation to fail with 64-bit R.
    //      
    // =========================
    // 1. Define all test data components as powers of a common base to make
    // it easy to accurately calculate the expected sum without adding any small
    // numbers together
    double start_num = 1125899906842624.0; // 2^50 = 1125899906842624
    double element_value = 0.03125; // 2^-5 = 0.03125
    double num_elements = 4096; // 2^12 = 4096
    double expected_total = 1125899906842752.0; // 2^50 + 2^7
    
    std::vector<double> input(num_elements, element_value);
    
    // Uncomment me if debugging test failure
    Rcerr << std::fixed << "Num elements: " << num_elements << "; Element value: " 
          << element_value << "; Exp. total: " << expected_total 
          << "; Sizeof(double): " << sizeof(double) << "\n";
    
    
    // Define acceptable tolerance
    double tolerance = element_value;

    // Do naiive and Kahan summation
    double naiive_total = start_num;
    double kahan_total = start_num;
    double kahan_compensation = 0.0;
    for(auto const& value: input) {
      naiive_total += value;
      add_element_kahan(kahan_total, value, kahan_compensation);
    }
    
    double naiive_diff = (expected_total - naiive_total);
    double kahan_diff = (expected_total - kahan_total);
    
    // Uncomment me if debugging test failure
    Rcerr << std::fixed << "Expected: " << expected_total << "; " << "Naiive: " 
          << naiive_total
          << "; Kahan: " << kahan_total << "; Naiive diff: " << naiive_diff
          << "; Kahan diff: " << kahan_diff << "; Tolerance: " << tolerance 
          << "\n";
    
    
    // Check naiive summation has unacceptable error
    expect_false(within_tolerance(naiive_total, expected_total, tolerance));
    // Check Kahan summation has acceptable error
    expect_true(within_tolerance(kahan_total, expected_total, tolerance));
    
  }
  
}
