// Enable C++11
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <math.h>

using namespace Rcpp;

//compute segment
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

  if (val1_start > val2_start) {
    if (val1_end >= val2_end) {
      // They are in the same order no bowtie
      // seg1 is above seg2
      // triangle of seg1
      topTriangle = 0.5*length*(val1_end-val1_start);
      // rectangle between seg1 and seg2
      topRectangle = length*(val1_start-val2_start);
      // triangle of seg2 (to be removed)
      bottomTriangle = 0.5*length*(val2_end-val2_start);
      return topTriangle+topRectangle-bottomTriangle;
    }
    else {
      //bowtie
      // lets make this really simple as the compiler
      // will combine the expressions as needed
      midPoint = (val1_start - val2_start) /
        ((val2_end - val2_start) - (val1_end - val1_start));

      midValue = val1_start + midPoint * (val1_end - val1_start);

      midPoint = midPoint * length;

      topTriangle = 0.5 * midPoint * (midValue - val1_start);
      topRectangle = midPoint * (val1_start - val2_start);
      bottomTriangle = 0.5 * midPoint * (midValue - val2_start);

      res = topTriangle + topRectangle - bottomTriangle;

      topTriangle = 0.5 * (length - midPoint) * (val2_end - midValue);
      topRectangle = 0; // midPoint*(val1_start-val2_start);
      bottomTriangle = 0.5 * (length - midPoint) * (val1_end - midValue);
      res += topTriangle + topRectangle - bottomTriangle;
      return res;
    }
  }
  else {
    if (val1_end > val2_end) {
      //bowtie
      // Find the point where they cross.
      // (Solution of linear equations)
      midPoint = (val2_start - val1_start) /
        ((val1_end - val1_start) - (val2_end - val2_start));

      midValue = val2_start + midPoint * (val2_end - val2_start);
      midPoint = midPoint * length;

      topTriangle = 0.5 * midPoint * (midValue - val2_start);
      topRectangle = midPoint * (val2_start - val1_start);
      bottomTriangle = 0.5 * midPoint * (midValue - val1_start);

      res = topTriangle + topRectangle - bottomTriangle;

      topTriangle = 0.5 * (length - midPoint) * (val1_end - midValue);
      topRectangle = 0; // midPoint*(val1_start-val2_start);
      bottomTriangle = 0.5 * (length - midPoint) * (val2_end - midValue);
      res += topTriangle + topRectangle - bottomTriangle;
      return res;
    }
    else { // same order
      // seg2 is above seg1
      // Triangle seg2 above seg1
      topTriangle = 0.5 * length * (val2_end - val2_start);
      // rectangle between seg2 and seg1
      topRectangle = length * (val2_start - val1_start);
      // Seg1 triangle to be removed
      bottomTriangle = 0.5 * length * (val1_end - val1_start);
      return topTriangle + topRectangle - bottomTriangle;
    }
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
  int index1, index2;

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

  double minLoc = std::min(loc1[0], loc2[0]);

  // warning area before loc2[0] is not well tested
  // As values outside of the range appear to be zero

  // Dealing with segments which are to the left of the region covered by both
  // segs
  if (loc2[0] < loc1[0]) {
    // loc2 starts before loc1 so lets deal with those segments first
    // Fix the position of Seg1 and then interate over Seg2 until we have all
    // of the segments of Seg2 before Seg1 starts.
    curSeg1Loc1 = minLoc;
    curSeg1Loc2 = minLoc;
    curSeg1Loc3 = loc1[0];
    curSeg1Val1 = 0;
    curSeg1Val2 = 0;

    // Set this value so we can update in the loop
    curSeg2Val2 = 0;
    for (index2 = 0; index2 < loc2.size(); index2++) {
      curSeg2Loc1 = loc2[index2];
      curSeg2Loc2 = loc2[index2] + binWidth2;
      if (index2 == loc2.size() - 1) {
        curSeg2Loc3 = maxLoc;
      }
      else {
        curSeg2Loc3 = loc2[index2 + 1];
      }
      curSeg2Val1 = curSeg2Val2;
      curSeg2Val2 = val2[index2];
      res += get_double_segment_constrained(
        curSeg1Loc1, curSeg1Loc2, curSeg1Loc3, curSeg1Val1, curSeg1Val2,
        curSeg2Loc1, curSeg2Loc2, curSeg2Loc3, curSeg2Val1, curSeg2Val2);

      if (curSeg1Loc1 > curSeg2Loc3) {
        break;
      }
    }
  }
  else
  {
    // loc2 starts before loc1 so lets deal with those segments first
    // Fix the position of Seg2 and then interate over Seg1 until we have all
    // of the segments of Seg1 before Seg2 starts.
    curSeg2Loc1 = minLoc;
    curSeg2Loc2 = minLoc;
    curSeg2Loc3 = loc2[0];
    curSeg2Val1 = 0;
    curSeg2Val2 = 0;

    // Set this value so we can update in the lopp
    curSeg1Val2 = 0;
    for (index1 = 0; index1 < loc1.size(); index1++) {
      curSeg1Loc1 = loc1[index1];
      curSeg1Loc2 = loc1[index1] + binWidth1;
      if (index1 == loc1.size() - 1) {
        curSeg1Loc3 = maxLoc;
      }
      else {
        curSeg1Loc3 = loc1[index1 + 1];
      }
      curSeg1Val1 = curSeg1Val2;
      curSeg1Val2 = val1[index1];
      res += get_double_segment_constrained(
        curSeg1Loc1, curSeg1Loc2, curSeg1Loc3, curSeg1Val1, curSeg1Val2,
        curSeg2Loc1, curSeg2Loc2, curSeg2Loc3, curSeg2Val1, curSeg2Val2);
      if (curSeg2Loc1 > curSeg1Loc3) {
        break;
      }
    }
  }

  // Add both the  overlapping sections and the non overlapping section on the right
  // Note we reiterate over the first few sections loc1
  // Could store where we are upto from above to save time
  // Reset Val counter
  curSeg1Val2 = 0;
  for (index1 = 0; index1 < loc1.size(); index1++) {
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
    for (index2 = secondStart; index2 < loc2.size(); index2++) {
      // Construct 3 sections for second seg.
      curSeg2Loc1=loc2[index2];
      curSeg2Loc2=loc2[index2] + binWidth2;
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
