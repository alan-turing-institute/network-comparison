// Enable C++11
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <math.h>

using namespace Rcpp;

inline double bowtie_area(double length, double val1_start, double val1_end,
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
inline double get_segment(double start, double end, double val1_start,
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
inline double get_segment_constrained(double seg1L1, double seg1L2,
                                      double seg2L1, double seg2L2,
                                      double seg1V1, double seg1V2,
                                      double seg2V1, double seg2V2)
{
  //We have a valid range
  double valStart1, valEnd1, valStart2, valEnd2;
  double start,end;
  start = std::max(seg1L1,seg2L1);
  end   = std::min(seg1L2,seg2L2);
  if (start < end) {
    valStart1 = seg1V1 + (seg1V2 - seg1V1)*(start - seg1L1)/(seg1L2 - seg1L1);
    valEnd1   = seg1V1 + (seg1V2 - seg1V1)*(end   - seg1L1)/(seg1L2 - seg1L1);
    valStart2 = seg2V1 + (seg2V2 - seg2V1)*(start - seg2L1)/(seg2L2 - seg2L1);
    valEnd2   = seg2V1 + (seg2V2 - seg2V1)*(end   - seg2L1)/(seg2L2 - seg2L1);

    return get_segment(start, end, valStart1, valEnd1, valStart2, valEnd2);
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


// Dealing with segments which are to the left of the region covered by both
// segs
inline double leftmost_segments(const NumericVector& loc1,
                                const NumericVector& loc2,
                                const NumericVector& val1,
                                const NumericVector& val2,
                                double binWidth1,
                                double maxLoc)
{
  double res = 0.0;

  assert(loc1[0] < loc2[0]);

  // Fix the position of Seg2 and then interate over Seg1 until we have all
  // of the segments of Seg1 before Seg2 starts.

  // are these all used?
  
  double curSeg2Loc1 = loc1[0];
  double curSeg2Loc2 = loc1[0];
  double curSeg2Loc3 = loc2[0];
  double curSeg2Val1 = 0;
  double curSeg2Val2 = 0;

  // Set this value so we can update in the lopp
  double curSeg1Val2 = 0;
  for (int index = 0; index < loc1.size(); index++) {
    double curSeg1Loc1 = loc1[index];
    double curSeg1Loc2 = loc1[index] + binWidth1;
    double curSeg1Loc3;
    if (index == loc1.size() - 1) {
      curSeg1Loc3 = maxLoc;
    }
    else {
      curSeg1Loc3 = loc1[index + 1];
    }
    double curSeg1Val1 = curSeg1Val2;
    curSeg1Val2 = val1[index];
    res += get_double_segment_constrained(
      curSeg1Loc1, curSeg1Loc2, curSeg1Loc3, curSeg1Val1, curSeg1Val2,
      curSeg2Loc1, curSeg2Loc2, curSeg2Loc3, curSeg2Val1, curSeg2Val2);

    if (curSeg2Loc1 > curSeg1Loc3) {
      break;
    }
  }

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
double NetEmdSmoothV2(NumericVector loc1, NumericVector val1, double binWidth1,
                      NumericVector loc2, NumericVector val2, double binWidth2)
{
  double res = 0;

  // Hist 1
  double curSeg1Loc1; // Start of the gradient section in Seg1
  double curSeg1Loc2; // End of the gradient section in Seg1
  double curSeg1Loc3; // End of the flat section in Seg1
  double curSeg1Val1; // Start value in Seg1
  double curSeg1Val2; // End value in Seg1

  // Hist 2
  double curSeg2Loc1; // Start of the gradient section in Seg2
  double curSeg2Loc2; // End of the gradient section in Seg2
  double curSeg2Loc3; // End of the flat section in Seg2
  double curSeg2Val1; // Start value in Seg2
  double curSeg2Val2; // End value in Seg2


  // Starting index for the second histogram
  double secondStart = 0;

  // Smallest and largest location values
  double maxLoc = std::max(loc1[loc1.size()-1] + binWidth1,
                           loc2[loc2.size()-1] + binWidth2);


  // warning area before loc2[0] is not well tested
  // As values outside of the range appear to be zero

  if (loc2[0] < loc1[0]) {
    res += leftmost_segments(loc2, loc1, val2, val1, binWidth2, maxLoc);
  }
  else
  {
    res += leftmost_segments(loc1, loc2, val1, val2, binWidth1, maxLoc);
  }

  // Add both the  overlapping sections and the non overlapping section on the right
  // Note we reiterate over the first few sections loc1
  // Could store where we are upto from above to save time
  // Reset Val counter
  curSeg1Val2 = 0;
  for (int index1 = 0; index1 < loc1.size(); index1++) {
    //  Get the three relevant locations
    // Start; end of linear section; end of flat section
    curSeg1Loc1 = loc1[index1];
    curSeg1Loc2 = loc1[index1] + binWidth1;
    // could pull this check outside of the loop with final case not sure if worth it
    if (index1 == loc1.size() - 1) {
      curSeg1Loc3 = maxLoc;
    }
    else {
      curSeg1Loc3 = loc1[index1 + 1];
    }
    // Update value to the start and end of the current section
    curSeg1Val1 = curSeg1Val2;
    curSeg1Val2 = val1[index1];
    // Setting up the previous value for the next loop
    // Could replace this loop with a while
    // but if so would need to be careful about overlaps
    if (secondStart == 0) {
      curSeg2Val2 = 0;
    }
    else {
      curSeg2Val2 = val2[secondStart - 1];
    }
    for (int index2 = secondStart; index2 < loc2.size(); index2++) {
      // Construct 3 sections for second seg.
      curSeg2Loc1 = loc2[index2];
      curSeg2Loc2 = loc2[index2] + binWidth2;
      if (index2 == loc2.size() - 1) {
        curSeg2Loc3 = maxLoc;
      }
      else {
        curSeg2Loc3 = loc2[index2 + 1];
      }
      //update values
      curSeg2Val1 = curSeg2Val2;
      curSeg2Val2 = val2[index2];
      // If this section is behind Seg1
      // Do not consider again
      if (curSeg2Loc3 < curSeg1Loc1) {
        secondStart = index2 + 1;
        continue;
      }
      // If current Seg2 is beyond Seg1 break out of loop
      res += get_double_segment_constrained(
        curSeg1Loc1, curSeg1Loc2, curSeg1Loc3, curSeg1Val1, curSeg1Val2,
        curSeg2Loc1, curSeg2Loc2, curSeg2Loc3, curSeg2Val1, curSeg2Val2);

      if (curSeg2Loc3 > curSeg1Loc3) {
        break;
      }
    }
  }
  return res;
}
