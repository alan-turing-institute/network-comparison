// Enable C++11
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <math.h>

using namespace Rcpp;

inline double bowtie_area(double length, double val1_start, double val1_end,
                          double val2_start, double val2_end);
inline double get_segment(double start, double end, double val1_start,
                          double val1_end, double val2_start, double val2_end);
inline double get_segment_constrained(double seg1L1, double seg1L2,
                                      double seg2L1, double seg2L2,
                                      double seg1V1, double seg1V2,
                                      double seg2V1, double seg2V2);

double get_double_segment_constrained(
  double seg1Loc1, double seg1Loc2, double seg1Loc3,
  double seg1Val1, double seg1Val2,
  double seg2Loc1, double seg2Loc2, double seg2Loc3,
  double seg2Val1, double seg2Val2);

inline double leftmost_segments(const NumericVector& loc1,
                                const NumericVector& loc2,
                                const NumericVector& val1,
                                const NumericVector& val2,
                                double binWidth1,
                                double maxLoc);

double NetEmdSmoothV2(NumericVector loc1, NumericVector val1, double binWidth1,
                      NumericVector loc2, NumericVector val2, double binWidth2);
