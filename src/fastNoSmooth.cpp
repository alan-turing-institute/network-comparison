// Enable C++11
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]


//' @title
//' Numerically safe addition using Kahan summation
////'
////' @param &sum Current accumulated sum. Updated by function.
////' @param &element Element to add to the accumulated sum
////' @param &compensation Current adjustment to compensate for floating point
////' summation error. Updated by function.
//'
//' @export
// [[Rcpp::export]]
void addElementKahan(double &sum, double element, double &compensation)
{
  // Uses Kahan summation algorithm to limit numerical error caused by adding
  // lots of very small things to a big thing to O(1) rather than O(N)
  // See: https://en.wikipedia.org/wiki/Kahan_summation_algorithm
  // Our implementation is the same as the Boost library implementation, with
  // the exception of any preventative measures taken to avoid unsafe math 
  // optimisations. The key compiler flag appears to be -ffast-math, which
  // allows the compiler to make optimisations that are not compliant with the
  // IEEE or ISO rules fdor math functions. We think that this must be set
  // esxplicitly, but in any case, we attempt to ensure it is uset by adding the
  // -fno-fast-math flag to PKG_CPPFLAGS in the Makevars file for the package.
  double compensatedElement = (element - compensation);
  double compensatedSum = sum + compensatedElement;
  compensation = ((compensatedSum - sum) - compensatedElement);
  sum = compensatedSum;
}

//' @title
//' Compute Earth Mover's Distance (EMD) between two Empirical Cumulative 
//' Density Functions (ECDFs)
////'
////' @param locations1 Locations for ECDF 1
////' @param values1 Cumulative masses for ECDF 1
////' @param locations2 Locations for ECDF 2
////' @param values2 Cumulative masses for ECDF 2
//'
//' @export
// [[Rcpp::export]]
double NetEmdConstant(NumericVector locations1, NumericVector values1, 
                      NumericVector locations2, NumericVector values2)
{
  double segmentStartLocation;
  double segmentArea;
  // Set start location of sweep below the minimum location across both ECDFs so
  // that we accumulate the full area between the two ECDFs
  if (locations1[0] < locations2[0])
  {
    segmentStartLocation = locations1[0] - 1.0;
  }
  else
  {
    segmentStartLocation = locations2[0] - 1.0;
  }
  double segmentValue1 = 0;
  double segmentValue2 = 0;
  int locationIndex1 = 0;
  int locationIndex2 = 0;
  double emd = 0;
  double maxValEcdf = 1.0;
  double compensation = 0.0;
  // We scan across all locations in both ECDFs, calculting the area of
  // each rectangular segment between the two ECDFs.
  // TODO: be worried about adding lots of small numbers
  while (true)
  {
    // Calculate the area of the next rectangular segment...
    if (locations1[locationIndex1] < locations2[locationIndex2])
    {
      // ...when next location is in ECDF 1
      segmentArea = (locations1[locationIndex1] - segmentStartLocation)
                    * std::abs(segmentValue1 - segmentValue2);
      addElementKahan(emd, segmentArea, compensation);
      segmentValue1 = values1[locationIndex1];
      segmentStartLocation = locations1[locationIndex1];
      locationIndex1 += 1;
      // If we've reached the end of ECDF 1, "short-circuit" the calculation
      // by stepping through the remaining points of ECDF 1 and then return
      if (locationIndex1 == locations1.size())
      {
        // Skip initialisation of loop variable as it is already set to correct
        // starting value
        for(; locationIndex2 < locations2.size(); locationIndex2++)
        {
          // No abs() in segment area calculation as we know max value will 
          // always be >= segment value for ECDF we are stepping through
          segmentArea = (locations2[locationIndex2] - segmentStartLocation) 
                        * (maxValEcdf - segmentValue2);
          addElementKahan(emd, segmentArea, compensation);
          segmentValue2 = values2[locationIndex2];
          segmentStartLocation = locations2[locationIndex2];
        }
        return emd;
      }
    }
    else
    {
      // ...when next location is in ECDF 2
      segmentArea = (locations2[locationIndex2] - segmentStartLocation) 
                    * std::abs(segmentValue1 - segmentValue2);
      addElementKahan(emd, segmentArea, compensation);
      segmentValue2 = values2[locationIndex2];
      segmentStartLocation = locations2[locationIndex2];
      locationIndex2 += 1;
      // If we've reached the end of ECDF 2, "short-circuit" the calculation
      // by stepping through the remaining points of ECDF 2 and then return
      if (locationIndex2 == locations2.size())
      {
        // Skip initialisation of loop variable as it is already set to correct
        // starting value
        for(; locationIndex1 < locations1.size(); locationIndex1++)
        {
          // No abs() in segment area calculation as we know max value will 
          // always be >= segment value for ECDF we are stepping through
          segmentArea = (locations1[locationIndex1] - segmentStartLocation)
                        * (maxValEcdf - segmentValue1);
          addElementKahan(emd, segmentArea, compensation);
          segmentValue1 = values1[locationIndex1];
          segmentStartLocation = locations1[locationIndex1];
        }
        return emd;
      }
    }
  }
}
