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
      emd += segmentArea;
      segmentValue1 = values1[locationIndex1];
      segmentStartLocation = locations1[locationIndex1];
      locationIndex1++;
      // If we've reached the end of ECDF 1, "short-circuit" the calculation
      // by stepping through the remaining points of ECDF 1 and then return
      if (locationIndex1 == locations1.size())
      {
        // Special first case as initial segmentStartLocation could be from
        // either ECDF. This lets us optimise the rest of the "short-circuit"
        // loop.
        segmentArea = (locations2[locationIndex2] - segmentStartLocation) 
                      * (maxValEcdf - segmentValue2);
        emd += segmentArea;
        locationIndex2++;
        // Skip initialisation of loop variable as it is already set to correct
        // starting value
        for(; locationIndex2 < locations2.size(); locationIndex2++)
        {
          // No abs() in segment area calculation as we know max value will 
          // always be >= segment value for ECDF we are stepping through
          // We stop using segmentStartLocation as, after initial case, this 
          // will always be the previous location in ECDF 2
          segmentArea = (locations2[locationIndex2] - locations2[locationIndex2-1]) 
                        * (maxValEcdf - values2[locationIndex2-1]);
          emd += segmentArea;
        }
        return emd;
      }
    }
    else
    {
      // ...when next location is in ECDF 2
      segmentArea = (locations2[locationIndex2] - segmentStartLocation) 
                    * std::abs(segmentValue1 - segmentValue2);
      emd += segmentArea;
      segmentValue2 = values2[locationIndex2];
      segmentStartLocation = locations2[locationIndex2];
      locationIndex2++;
      // If we've reached the end of ECDF 2, "short-circuit" the calculation
      // by stepping through the remaining points of ECDF 2 and then return
      if (locationIndex2 == locations2.size())
      {
        // Special first case as initial segmentStartLocation could be from
        // either ECDF. This lets us optimise the rest of the "short-circuit"
        // loop.
        segmentArea = (locations1[locationIndex1] - segmentStartLocation)
                      * (maxValEcdf - segmentValue1);
        emd += segmentArea;
        locationIndex1++;
        // Skip initialisation of loop variable as it is already set to correct
        // starting value
        for(; locationIndex1 < locations1.size(); locationIndex1++)
        {
          // No abs() in segment area calculation as we know max value will 
          // always be >= segment value for ECDF we are stepping through
          // We stop using segmentStartLocation as, after initial case, this 
          // will always be the previous location in ECDF 1
          segmentArea = (locations1[locationIndex1] - locations1[locationIndex1-1])
                        * (maxValEcdf - values1[locationIndex1-1]);
          emd += segmentArea;
        }
        return emd;
      }
    }
  }
}
