// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// logistic_regression
NumericVector logistic_regression(NumericMatrix explain, NumericVector response);
RcppExport SEXP _StatComp22011_logistic_regression(SEXP explainSEXP, SEXP responseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type explain(explainSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type response(responseSEXP);
    rcpp_result_gen = Rcpp::wrap(logistic_regression(explain, response));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_StatComp22011_logistic_regression", (DL_FUNC) &_StatComp22011_logistic_regression, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_StatComp22011(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
