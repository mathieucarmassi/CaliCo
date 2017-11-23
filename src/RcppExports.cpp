// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// MetropolisHastingsCpp
List MetropolisHastingsCpp(Function model, int Ngibbs, int Nmh, arma::vec theta_init, arma::vec k, arma::mat SIGMA, arma::vec Yf, arma::vec binf, arma::vec bsup, Function LogTest, int stream);
RcppExport SEXP _calibrationCode_MetropolisHastingsCpp(SEXP modelSEXP, SEXP NgibbsSEXP, SEXP NmhSEXP, SEXP theta_initSEXP, SEXP kSEXP, SEXP SIGMASEXP, SEXP YfSEXP, SEXP binfSEXP, SEXP bsupSEXP, SEXP LogTestSEXP, SEXP streamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Function >::type model(modelSEXP);
    Rcpp::traits::input_parameter< int >::type Ngibbs(NgibbsSEXP);
    Rcpp::traits::input_parameter< int >::type Nmh(NmhSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta_init(theta_initSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type SIGMA(SIGMASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Yf(YfSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type binf(binfSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type bsup(bsupSEXP);
    Rcpp::traits::input_parameter< Function >::type LogTest(LogTestSEXP);
    Rcpp::traits::input_parameter< int >::type stream(streamSEXP);
    rcpp_result_gen = Rcpp::wrap(MetropolisHastingsCpp(model, Ngibbs, Nmh, theta_init, k, SIGMA, Yf, binf, bsup, LogTest, stream));
    return rcpp_result_gen;
END_RCPP
}
// invMat
arma::mat invMat(arma::mat V);
RcppExport SEXP _calibrationCode_invMat(SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(invMat(V));
    return rcpp_result_gen;
END_RCPP
}
// FlushCPP
void FlushCPP();
RcppExport SEXP _calibrationCode_FlushCPP() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    FlushCPP();
    return R_NilValue;
END_RCPP
}
// MetropolisHastingsCppD
List MetropolisHastingsCppD(Function model, int Ngibbs, int Nmh, arma::vec theta_init, arma::vec k, arma::mat SIGMA, arma::vec Yf, arma::vec binf, arma::vec bsup, Function LogTest);
RcppExport SEXP _calibrationCode_MetropolisHastingsCppD(SEXP modelSEXP, SEXP NgibbsSEXP, SEXP NmhSEXP, SEXP theta_initSEXP, SEXP kSEXP, SEXP SIGMASEXP, SEXP YfSEXP, SEXP binfSEXP, SEXP bsupSEXP, SEXP LogTestSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Function >::type model(modelSEXP);
    Rcpp::traits::input_parameter< int >::type Ngibbs(NgibbsSEXP);
    Rcpp::traits::input_parameter< int >::type Nmh(NmhSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta_init(theta_initSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type SIGMA(SIGMASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Yf(YfSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type binf(binfSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type bsup(bsupSEXP);
    Rcpp::traits::input_parameter< Function >::type LogTest(LogTestSEXP);
    rcpp_result_gen = Rcpp::wrap(MetropolisHastingsCppD(model, Ngibbs, Nmh, theta_init, k, SIGMA, Yf, binf, bsup, LogTest));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_calibrationCode_MetropolisHastingsCpp", (DL_FUNC) &_calibrationCode_MetropolisHastingsCpp, 11},
    {"_calibrationCode_invMat", (DL_FUNC) &_calibrationCode_invMat, 1},
    {"_calibrationCode_FlushCPP", (DL_FUNC) &_calibrationCode_FlushCPP, 0},
    {"_calibrationCode_MetropolisHastingsCppD", (DL_FUNC) &_calibrationCode_MetropolisHastingsCppD, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_calibrationCode(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
