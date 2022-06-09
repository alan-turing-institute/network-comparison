// Enable C++11
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <math.h>
#include "emd_fast_no_smoothing.h" // add_element_kahan()
#include "fastSmoothV2.h"

using namespace Rcpp;

double bowtie_area(double length, double val1_start, double val1_end,
                   double val2_start, double val2_end)
{
  double midPoint = (val1_start - val2_start) /
    ((val2_end - val2_start) - (val1_end - val1_start));

  const double midValue = val1_start + midPoint * (val1_end - val1_start);

  midPoint = midPoint * length;

  double topTriangle = 0.5 * midPoint * (midValue - val1_start);
  double topRectangle = midPoint * (val1_start - val2_start);
  double bottomTriangle = 0.5 * midPoint * (midValue - val2_start);

  double res = topTriangle + topRectangle - bottomTriangle;

  topTriangle = 0.5 * (length - midPoint) * (val2_end - midValue);
  topRectangle = 0; // midPoint*(val1_start-val2_start);
  bottomTriangle = 0.5 * (length - midPoint) * (val1_end - midValue);

  res += topTriangle + topRectangle - bottomTriangle;
  return res;
}

// Compute the unsigned area between two line segments
// assumes that val1_end > val1_start and val2_end > val2_start
double get_segment(double start, double end, double val1_start,
                   double val1_end, double val2_start, double val2_end)
{
  const double length = end - start;

  double topTriangle;
  double topRectangle;
  double bottomTriangle;
  double midPoint;
  double midValue;
  double res = 0;

  bool both_differences_positive = val1_start > val2_start && val1_end >= val2_end;
  bool both_differences_negative = val1_start <= val2_start && val1_end <= val2_end;

  if (both_differences_positive || both_differences_negative)
  {
    // They are in the same order: no bowtie
    // triangle of seg1
    topTriangle = 0.5 * length * (val1_end - val1_start);
    // rectangle between seg1 and seg2
    topRectangle = length * (val1_start - val2_start);
    // triangle of seg2 (to be removed)
    bottomTriangle = 0.5 * length * (val2_end - val2_start);

    const double sign = both_differences_positive?1.0:-1.0;
    return sign * (topTriangle + topRectangle - bottomTriangle);
  }
  else if (val1_start > val2_start) { // bowtie, first case
    return bowtie_area(length, val1_start, val1_end, val2_start, val2_end);
  }
  else { // bowtie, second case
    return bowtie_area(length, val2_start, val2_end, val1_start, val1_end);
  }
}

// cut down and compute segment
double get_segment_constrained(double seg1L1, double seg1L2,
                               double seg2L1, double seg2L2,
                               double seg1V1, double seg1V2,
                               double seg2V1, double seg2V2)
{
  //We have a valid range
  double valStart1, valEnd1, valStart2, valEnd2;
  double start,end;
  double result;
  start = std::max(seg1L1,seg2L1);
  end   = std::min(seg1L2,seg2L2);
  if (start < end) {
    valStart1 = seg1V1 + (seg1V2 - seg1V1)*(start - seg1L1)/(seg1L2 - seg1L1);
    valEnd1   = seg1V1 + (seg1V2 - seg1V1)*(end   - seg1L1)/(seg1L2 - seg1L1);
    valStart2 = seg2V1 + (seg2V2 - seg2V1)*(start - seg2L1)/(seg2L2 - seg2L1);
    valEnd2   = seg2V1 + (seg2V2 - seg2V1)*(end   - seg2L1)/(seg2L2 - seg2L1);
    result = get_segment(start, end, valStart1, valEnd1, valStart2, valEnd2);
    return result;
  }
  else {
    return 0;
  }
}

double get_double_segment_constrained(
  double seg1Loc1, double seg1Loc2, double seg1Loc3,
  double seg1Val1, double seg1Val2,
  double seg2Loc1, double seg2Loc2, double seg2Loc3,
  double seg2Val1, double seg2Val2)
{
  double res = 0;

  // compare the linear section with the linear section
  res += get_segment_constrained(seg1Loc1, seg1Loc2, seg2Loc1, seg2Loc2,
                                 seg1Val1, seg1Val2, seg2Val1, seg2Val2);

  // compare the linear section with the flat section
  // This could be easily special cased (saving ~1 if statements ).
  res += get_segment_constrained(seg1Loc1, seg1Loc2, seg2Loc2, seg2Loc3,
                                 seg1Val1, seg1Val2, seg2Val2, seg2Val2);

  // compare the flat section with the linear section
  // This could be easily special cased (saving ~1 if statements ).
  res += get_segment_constrained(seg1Loc2, seg1Loc3, seg2Loc1, seg2Loc2,
                                 seg1Val2, seg1Val2, seg2Val1, seg2Val2);

  // compare the flat section with the flat section
  // This could be easily special cased (saving ~2 if statements ).
  res += get_segment_constrained(seg1Loc2, seg1Loc3, seg2Loc2, seg2Loc3,
                                 seg1Val2, seg1Val2, seg2Val2, seg2Val2);

  return res;
}

//' @title
//' Compute EMD
////'
////' @param loc1 numeric vector.
////' @param val1 numeric vector.
////' @param loc2 numeric vector.
////' @param val2 numeric vector.
//'
//' @export
// [[Rcpp::export]]
double netemd_smooth(NumericVector loc1, NumericVector val1, double binWidth1,
                      NumericVector loc2, NumericVector val2, double binWidth2)
{
  OverlappingSegments segs(loc1, loc2, binWidth1, binWidth2);

  double res = 0.0;
  double compensation = 0.0;

  for (OverlappingSegments::iterator it = segs.begin(), endIt = segs.end();
       it != endIt; ++it)
  {
    // The OverlappingSegments iterator returns pairs of the left
    // indices of overlapping endpoints: it->first is the index into
    // loc1, it->second is the index into loc2.  When the smallest
    // element in one of these vectors is larger than the current
    // element of the other, an 'index' of -1 is returned.

    // Hist 1
    // Start of the gradient section in Seg1
    double curSeg1Loc1 = segs.loc1_left(it->first);

    // End of the gradient section in Seg1
    double curSeg1Loc2 = segs.loc1_mid(it->first);

    // End of the flat section in Seg1
    double curSeg1Loc3 = segs.loc1_right(it->first);

    // Start and end values in Seg1: val1 gives the values at *right*
    // endpoints of the segments.  A value of 0.0 is used before the
    // first segment.
    double curSeg1Val1 = (it->first > 0) ? val1[it->first - 1] : 0.0;
    double curSeg1Val2 = (it->first >= 0) ? val1[it->first] : 0.0;

    // Hist 2
    // Start of the gradient section in Seg1
    double curSeg2Loc1 = segs.loc2_left(it->second);

    // End of the gradient section in Seg1
    double curSeg2Loc2 = segs.loc2_mid(it->second);

    // End of the flat section in Seg1
    double curSeg2Loc3 = segs.loc2_right(it->second);

    // Start and end values in Seg2: val2 gives the values at *right*
    // endpoints of the segments.  A value of 0.0 is used before the
    // first segment.
    double curSeg2Val1 = (it->second > 0) ? val2[it->second - 1] : 0.0;
    double curSeg2Val2 = (it->second >= 0) ? val2[it->second] : 0.0;

    double element = get_double_segment_constrained(
      curSeg1Loc1, curSeg1Loc2, curSeg1Loc3, curSeg1Val1, curSeg1Val2,
      curSeg2Loc1, curSeg2Loc2, curSeg2Loc3, curSeg2Val1, curSeg2Val2);

    add_element_kahan(res, element, compensation);
  }

  return res;
}
