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
//' Compute Earth Mover's Distance (EMD) between two Empirical Cumulative Mass
//' Functions (ECMFs)
////'
////' @param locations1 Locations for ECMF 1
////' @param values1 Cumulative masses for ECMF 1
////' @param locations2 Locations for ECMF 2
////' @param values2 Cumulative masses for ECMF 2
//'
//' @export
// [[Rcpp::export]]
double NetEmdConstant(NumericVector locations1, NumericVector values1, 
                      NumericVector locations2, NumericVector values2)
{
  double currentLocation;
  double segmentArea;
  // Set start location of sweep below the minimum location across both ECMFs so
  // that we accumulate the full area between the two ECMFs
  if (locations1[0] < locations2[0])
  {
    currentLocation = locations1[0] - 1.0;
  }
  else
  {
    currentLocation = locations2[0] - 1.0;
  }
  double currentValue1 = 0;
  double currentValue2 = 0;
  int locationIndex1 = 0;
  int locationIndex2 = 0;
  double emd = 0;
  // We scan across all locations in both ECMFs, calculting the area of
  // each rectangular segment between the two ECMFs.
  // NOTE: We stop when we hit the last location on either ECMF and handle the
  // "tail" beyond this as a special case
  // TODO: be worried about adding lots of small numbers
  while (locationIndex1 < locations1.size() ||
         locationIndex2 < locations2.size())
  {
    // Stop if we reach the end of either ECMF
    if (locationIndex1 == locations1.size())
    {
      // We are beyond the last location of ECMF 1
      emd+=(locations2[locationIndex2]-currentLocation)*std::abs(currentValue1-currentValue2);
      currentValue2=values2[locationIndex2];
      currentLocation=locations2[locationIndex2];
      locationIndex2 += 1;
    }
    else if (locationIndex2 == locations2.size())
    {
      // We are beyond the last location of ECMF 2
      emd+=(locations1[locationIndex1]-currentLocation)*std::abs(currentValue2-currentValue1);
      currentValue1=values1[locationIndex1];
      currentLocation=locations1[locationIndex1];
      locationIndex1 += 1;
    }
    // Calculate the area of the next rectangular segment...
    else if (locations1[locationIndex1] < locations2[locationIndex2])
    {
      // ...when next location is in ECMF 1
      segmentArea = (locations1[locationIndex1] - currentLocation)
                    * std::abs(currentValue1 - currentValue2);
      emd += segmentArea;
      currentValue1 = values1[locationIndex1];
      currentLocation = locations1[locationIndex1];
      locationIndex1 += 1;
    }
    else
    {
      // ...when next location is in ECMF 2
      segmentArea = (locations2[locationIndex2] - currentLocation) 
                    * std::abs(currentValue1 - currentValue2);
      emd += segmentArea;
      currentValue2 = values2[locationIndex2];
      currentLocation = locations2[locationIndex2];
      locationIndex2 += 1;
    }
  }
  return emd;
}
