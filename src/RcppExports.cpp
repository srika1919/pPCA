// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// colMSD_dgc
NumericVector colMSD_dgc(S4 mat, NumericVector m);
RcppExport SEXP _pPCA_colMSD_dgc(SEXP matSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type mat(matSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(colMSD_dgc(mat, m));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_pPCA_colMSD_dgc", (DL_FUNC) &_pPCA_colMSD_dgc, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_pPCA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
