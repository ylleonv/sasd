// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _anRpackage_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_mimod();
RcppExport SEXP _rcpp_module_boot_NumEx();
RcppExport SEXP _rcpp_module_boot_yada();
RcppExport SEXP _rcpp_module_boot_stdVector();

static const R_CallMethodDef CallEntries[] = {
    {"_anRpackage_rcpp_hello_world", (DL_FUNC) &_anRpackage_rcpp_hello_world, 0},
    {"_rcpp_module_boot_mimod", (DL_FUNC) &_rcpp_module_boot_mimod, 0},
    {"_rcpp_module_boot_NumEx", (DL_FUNC) &_rcpp_module_boot_NumEx, 0},
    {"_rcpp_module_boot_yada", (DL_FUNC) &_rcpp_module_boot_yada, 0},
    {"_rcpp_module_boot_stdVector", (DL_FUNC) &_rcpp_module_boot_stdVector, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_anRpackage(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
