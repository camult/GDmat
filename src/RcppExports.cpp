// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// DRM
arma::mat DRM(arma::mat X);
RcppExport SEXP _GDmat_DRM(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(DRM(X));
    return rcpp_result_gen;
END_RCPP
}
// GRM
arma::mat GRM(arma::mat X, double smallValue);
RcppExport SEXP _GDmat_GRM(SEXP XSEXP, SEXP smallValueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type smallValue(smallValueSEXP);
    rcpp_result_gen = Rcpp::wrap(GRM(X, smallValue));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GDmat_DRM", (DL_FUNC) &_GDmat_DRM, 1},
    {"_GDmat_GRM", (DL_FUNC) &_GDmat_GRM, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_GDmat(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
