// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// counts_from_observations
NumericMatrix counts_from_observations(NumericMatrix features);
RcppExport SEXP _netdist_counts_from_observations(SEXP featuresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type features(featuresSEXP);
    rcpp_result_gen = Rcpp::wrap(counts_from_observations(features));
    return rcpp_result_gen;
END_RCPP
}
// emd_fast_no_smoothing
double emd_fast_no_smoothing(NumericVector locations1, NumericVector values1, NumericVector locations2, NumericVector values2);
RcppExport SEXP _netdist_emd_fast_no_smoothing(SEXP locations1SEXP, SEXP values1SEXP, SEXP locations2SEXP, SEXP values2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type locations1(locations1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type values1(values1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type locations2(locations2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type values2(values2SEXP);
    rcpp_result_gen = Rcpp::wrap(emd_fast_no_smoothing(locations1, values1, locations2, values2));
    return rcpp_result_gen;
END_RCPP
}
// NetEmdSmooth
double NetEmdSmooth(NumericVector loc1, NumericVector val1, double binWidth1, NumericVector loc2, NumericVector val2, double binWidth2);
RcppExport SEXP _netdist_NetEmdSmooth(SEXP loc1SEXP, SEXP val1SEXP, SEXP binWidth1SEXP, SEXP loc2SEXP, SEXP val2SEXP, SEXP binWidth2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type loc1(loc1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type val1(val1SEXP);
    Rcpp::traits::input_parameter< double >::type binWidth1(binWidth1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type loc2(loc2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type val2(val2SEXP);
    Rcpp::traits::input_parameter< double >::type binWidth2(binWidth2SEXP);
    rcpp_result_gen = Rcpp::wrap(NetEmdSmooth(loc1, val1, binWidth1, loc2, val2, binWidth2));
    return rcpp_result_gen;
END_RCPP
}
// NetEmdSmoothV2
double NetEmdSmoothV2(NumericVector loc1, NumericVector val1, double binWidth1, NumericVector loc2, NumericVector val2, double binWidth2);
RcppExport SEXP _netdist_NetEmdSmoothV2(SEXP loc1SEXP, SEXP val1SEXP, SEXP binWidth1SEXP, SEXP loc2SEXP, SEXP val2SEXP, SEXP binWidth2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type loc1(loc1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type val1(val1SEXP);
    Rcpp::traits::input_parameter< double >::type binWidth1(binWidth1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type loc2(loc2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type val2(val2SEXP);
    Rcpp::traits::input_parameter< double >::type binWidth2(binWidth2SEXP);
    rcpp_result_gen = Rcpp::wrap(NetEmdSmoothV2(loc1, val1, binWidth1, loc2, val2, binWidth2));
    return rcpp_result_gen;
END_RCPP
}
// NetEmdSmoothV2_old
double NetEmdSmoothV2_old(NumericVector loc1, NumericVector val1, double binWidth1, NumericVector loc2, NumericVector val2, double binWidth2);
RcppExport SEXP _netdist_NetEmdSmoothV2_old(SEXP loc1SEXP, SEXP val1SEXP, SEXP binWidth1SEXP, SEXP loc2SEXP, SEXP val2SEXP, SEXP binWidth2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type loc1(loc1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type val1(val1SEXP);
    Rcpp::traits::input_parameter< double >::type binWidth1(binWidth1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type loc2(loc2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type val2(val2SEXP);
    Rcpp::traits::input_parameter< double >::type binWidth2(binWidth2SEXP);
    rcpp_result_gen = Rcpp::wrap(NetEmdSmoothV2_old(loc1, val1, binWidth1, loc2, val2, binWidth2));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP run_testthat_tests();

static const R_CallMethodDef CallEntries[] = {
    {"_netdist_counts_from_observations", (DL_FUNC) &_netdist_counts_from_observations, 1},
    {"_netdist_emd_fast_no_smoothing", (DL_FUNC) &_netdist_emd_fast_no_smoothing, 4},
    {"_netdist_NetEmdSmooth", (DL_FUNC) &_netdist_NetEmdSmooth, 6},
    {"_netdist_NetEmdSmoothV2", (DL_FUNC) &_netdist_NetEmdSmoothV2, 6},
    {"_netdist_NetEmdSmoothV2_old", (DL_FUNC) &_netdist_NetEmdSmoothV2_old, 6},
    {"run_testthat_tests", (DL_FUNC) &run_testthat_tests, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_netdist(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
