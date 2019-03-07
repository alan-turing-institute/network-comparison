#include <Rcpp.h>
using namespace Rcpp;

void add_element_kahan(double &sum, double element, double &compensation);
double emd_fast_no_smoothing(NumericVector locations1, NumericVector values1, 
                             NumericVector locations2, NumericVector values2);