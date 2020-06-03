// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// simulate_gaussian
NumericVector simulate_gaussian(int t_len, int n_ts, int sim_runs, Rcpp::NumericVector gset, int n, Rcpp::NumericVector sigma_vec, int deriv_order, bool epidem);
RcppExport SEXP _multiscale_simulate_gaussian(SEXP t_lenSEXP, SEXP n_tsSEXP, SEXP sim_runsSEXP, SEXP gsetSEXP, SEXP nSEXP, SEXP sigma_vecSEXP, SEXP deriv_orderSEXP, SEXP epidemSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type t_len(t_lenSEXP);
    Rcpp::traits::input_parameter< int >::type n_ts(n_tsSEXP);
    Rcpp::traits::input_parameter< int >::type sim_runs(sim_runsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type gset(gsetSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma_vec(sigma_vecSEXP);
    Rcpp::traits::input_parameter< int >::type deriv_order(deriv_orderSEXP);
    Rcpp::traits::input_parameter< bool >::type epidem(epidemSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_gaussian(t_len, n_ts, sim_runs, gset, n, sigma_vec, deriv_order, epidem));
    return rcpp_result_gen;
END_RCPP
}
// compute_multiple_statistics
Rcpp::List compute_multiple_statistics(int t_len, int n_ts, Rcpp::NumericMatrix data, Rcpp::NumericVector gset, int n, Rcpp::NumericVector sigma_vec, int deriv_order, bool epidem);
RcppExport SEXP _multiscale_compute_multiple_statistics(SEXP t_lenSEXP, SEXP n_tsSEXP, SEXP dataSEXP, SEXP gsetSEXP, SEXP nSEXP, SEXP sigma_vecSEXP, SEXP deriv_orderSEXP, SEXP epidemSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type t_len(t_lenSEXP);
    Rcpp::traits::input_parameter< int >::type n_ts(n_tsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type gset(gsetSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma_vec(sigma_vecSEXP);
    Rcpp::traits::input_parameter< int >::type deriv_order(deriv_orderSEXP);
    Rcpp::traits::input_parameter< bool >::type epidem(epidemSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_multiple_statistics(t_len, n_ts, data, gset, n, sigma_vec, deriv_order, epidem));
    return rcpp_result_gen;
END_RCPP
}
// compute_statistics
Rcpp::List compute_statistics(int t_len, Rcpp::NumericVector data, Rcpp::NumericVector gset, int n, double sigma, int deriv_order);
RcppExport SEXP _multiscale_compute_statistics(SEXP t_lenSEXP, SEXP dataSEXP, SEXP gsetSEXP, SEXP nSEXP, SEXP sigmaSEXP, SEXP deriv_orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type t_len(t_lenSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type gset(gsetSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type deriv_order(deriv_orderSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_statistics(t_len, data, gset, n, sigma, deriv_order));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_multiscale_simulate_gaussian", (DL_FUNC) &_multiscale_simulate_gaussian, 8},
    {"_multiscale_compute_multiple_statistics", (DL_FUNC) &_multiscale_compute_multiple_statistics, 8},
    {"_multiscale_compute_statistics", (DL_FUNC) &_multiscale_compute_statistics, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_multiscale(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
