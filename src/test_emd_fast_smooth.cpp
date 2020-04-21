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
#include "fastSmoothV2.h"
#include "emd_fast_no_smoothing.h"
#include <math.h>
#include <vector>
#include <testthat.h>

// Helper function to test tolerance
bool within_toleranceV2(double actual, double expected, double tolerance) {
  if(actual > expected) {
    return ((actual - expected) <= tolerance);
  }
  else {
    return ((expected - actual) <= tolerance);
  }
}

double simpleSlowArea(double startx,double endx,double starty1,double endy1,double starty2,double endy2)
{
        // Making this step size smaller 
        double step = (endx-startx)/100000000.0;
        double curX;
        double curY1;
        double curY2;
        double res = 0;
        for (int i=0;i<100000000;i++)
        {
                curX = startx + i*step;
                curY1 = starty1 +(endy1-starty1)*i/100000000.0;
                curY2 = starty2 +(endy2-starty2)*i/100000000.0;
                res += step*std::abs(curY1-curY2);
        }
        return res;
}

void runSegmentConstraintTest(double start,double end,double val1_start,double val1_end,double val2_start,double val2_end)
{
    double tempVal1;
    double tempVal2;
    tempVal1 = get_segment_constrained(start,end,start,end,val1_start,val1_end,val2_start,val2_end);
    tempVal2 = simpleSlowArea(start,end, val1_start,val1_end,val2_start,val2_end);
    std::cout << "\n";
    std::cout << "segment constrained " << tempVal1 << " simpleResult " << tempVal2 << "\n";
    expect_true(within_toleranceV2(tempVal1,tempVal2,0.0001));
}


void runSegmentTest(double start,double end,double val1_start,double val1_end,double val2_start,double val2_end)
{
    double tempVal1;
    double tempVal2;
    tempVal1 = get_segment(start,end, val1_start,val1_end,val2_start,val2_end);
    tempVal2 = simpleSlowArea(start,end, val1_start,val1_end,val2_start,val2_end);
    std::cout << "\n";
    std::cout << "segment test " << tempVal1 << " simpleResult " << tempVal2 << "\n";
    expect_true(within_toleranceV2(tempVal1,tempVal2,0.0001));
}


context("emd_fast_smoothing segment constrain simple") {
  test_that("emd_fast_smoothing segment constrain simple") {
    // Two upward linear segments
    runSegmentConstraintTest(0.0,1.0,0.0,1.0,0.0,1.0);
    // One upward one down linear segments
    runSegmentConstraintTest(0.0,1.0,0.0,1.0,1.0,0.0);
    runSegmentConstraintTest(0.0,1.0,1.0,0.0,0.0,1.0);
    // Two down linear segments
    runSegmentConstraintTest(0.0,1.0,1.0,0.0,1.0,0.0);
    // One flat one up segments
    runSegmentConstraintTest(0.0,1.0,0.0,0.0,0.0,1.0);
    runSegmentConstraintTest(0.0,1.0,0.0,1.0,0.0,0.0);
    // One flat one down segments
    runSegmentConstraintTest(0.0,1.0,1.0,0.0,0.0,0.0);
    runSegmentConstraintTest(0.0,1.0,0.0,0.0,1.0,0.0);
    // Different gradients segments
    runSegmentConstraintTest(0.0,1.0,0.0,3.0,0.0,0.0);
    runSegmentConstraintTest(0.0,1.0,0.0,0.0,0.0,3.0);
    // Different gradients segments
    runSegmentConstraintTest(0.0,1.0,2.0,4.0,1.0,2.0);
    runSegmentConstraintTest(0.0,1.0,1.0,2.0,2.0,3.0);
}}

context("emd_fast_smoothing segment full") {
  test_that("emd_fast_smoothing segment full") {
    // Two upward linear segments
    runSegmentTest(0.0,1.0,0.0,1.0,0.0,1.0);
    // One upward one down linear segments
    runSegmentTest(0.0,1.0,0.0,1.0,1.0,0.0);
    runSegmentTest(0.0,1.0,1.0,0.0,0.0,1.0);
    // Two down linear segments
    runSegmentTest(0.0,1.0,1.0,0.0,1.0,0.0);
    // One flat one up segments
    runSegmentTest(0.0,1.0,0.0,0.0,0.0,1.0);
    runSegmentTest(0.0,1.0,0.0,1.0,0.0,0.0);
    // One flat one down segments
    runSegmentTest(0.0,1.0,1.0,0.0,0.0,0.0);
    runSegmentTest(0.0,1.0,0.0,0.0,1.0,0.0);
    // Different gradients segments
    runSegmentTest(0.0,1.0,0.0,3.0,0.0,0.0);
    runSegmentTest(0.0,1.0,0.0,0.0,0.0,3.0);
    // Different gradients segments
    runSegmentTest(0.0,1.0,2.0,4.0,1.0,2.0);
    runSegmentTest(0.0,1.0,1.0,2.0,2.0,3.0);
}}


