// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// covNonSep
Eigen::MatrixXd covNonSep(const Eigen::MatrixXd DS, const Eigen::MatrixXd DT, double c, double thetas, double thetat);
RcppExport SEXP _STplm_covNonSep(SEXP DSSEXP, SEXP DTSEXP, SEXP cSEXP, SEXP thetasSEXP, SEXP thetatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type DS(DSSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type DT(DTSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type thetas(thetasSEXP);
    Rcpp::traits::input_parameter< double >::type thetat(thetatSEXP);
    rcpp_result_gen = Rcpp::wrap(covNonSep(DS, DT, c, thetas, thetat));
    return rcpp_result_gen;
END_RCPP
}
// covNonSep2
Eigen::MatrixXd covNonSep2(const Eigen::MatrixXd DS, const Eigen::MatrixXd DT, double c, double thetat, double thetas);
RcppExport SEXP _STplm_covNonSep2(SEXP DSSEXP, SEXP DTSEXP, SEXP cSEXP, SEXP thetatSEXP, SEXP thetasSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type DS(DSSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type DT(DTSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type thetat(thetatSEXP);
    Rcpp::traits::input_parameter< double >::type thetas(thetasSEXP);
    rcpp_result_gen = Rcpp::wrap(covNonSep2(DS, DT, c, thetat, thetas));
    return rcpp_result_gen;
END_RCPP
}
// covNonSep3
Eigen::MatrixXd covNonSep3(const Eigen::MatrixXd DS, const Eigen::MatrixXd DT, double c, double thetat, double thetas, double theta3);
RcppExport SEXP _STplm_covNonSep3(SEXP DSSEXP, SEXP DTSEXP, SEXP cSEXP, SEXP thetatSEXP, SEXP thetasSEXP, SEXP theta3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type DS(DSSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type DT(DTSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type thetat(thetatSEXP);
    Rcpp::traits::input_parameter< double >::type thetas(thetasSEXP);
    Rcpp::traits::input_parameter< double >::type theta3(theta3SEXP);
    rcpp_result_gen = Rcpp::wrap(covNonSep3(DS, DT, c, thetat, thetas, theta3));
    return rcpp_result_gen;
END_RCPP
}
// covSep
Eigen::MatrixXd covSep(const Eigen::MatrixXd DS, const Eigen::MatrixXd DT, double c, double thetas, double thetat);
RcppExport SEXP _STplm_covSep(SEXP DSSEXP, SEXP DTSEXP, SEXP cSEXP, SEXP thetasSEXP, SEXP thetatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type DS(DSSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type DT(DTSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type thetas(thetasSEXP);
    Rcpp::traits::input_parameter< double >::type thetat(thetatSEXP);
    rcpp_result_gen = Rcpp::wrap(covSep(DS, DT, c, thetas, thetat));
    return rcpp_result_gen;
END_RCPP
}
// loglikeCpp
Eigen::MatrixXd loglikeCpp(const Eigen::VectorXd theta0, const int CovStructure, const Eigen::MatrixXd loctim, const Eigen::MatrixXd DS, const Eigen::MatrixXd DT, const Eigen::VectorXd Z, const Eigen::VectorXd truncate, const double lambda, const int pterm, const Eigen::VectorXd DDnew);
RcppExport SEXP _STplm_loglikeCpp(SEXP theta0SEXP, SEXP CovStructureSEXP, SEXP loctimSEXP, SEXP DSSEXP, SEXP DTSEXP, SEXP ZSEXP, SEXP truncateSEXP, SEXP lambdaSEXP, SEXP ptermSEXP, SEXP DDnewSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type theta0(theta0SEXP);
    Rcpp::traits::input_parameter< const int >::type CovStructure(CovStructureSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type loctim(loctimSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type DS(DSSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type DT(DTSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type truncate(truncateSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const int >::type pterm(ptermSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type DDnew(DDnewSEXP);
    rcpp_result_gen = Rcpp::wrap(loglikeCpp(theta0, CovStructure, loctim, DS, DT, Z, truncate, lambda, pterm, DDnew));
    return rcpp_result_gen;
END_RCPP
}
// logProfileCpp
Eigen::MatrixXd logProfileCpp(const Eigen::VectorXd theta0, const int CovStructure, const Eigen::MatrixXd loctim, const Eigen::MatrixXd DS, const Eigen::MatrixXd DT, const Eigen::MatrixXd Xprime, const Eigen::VectorXd Zprime, const Eigen::VectorXd truncate, const double lambda, const int pterm, const Eigen::VectorXd DDnew);
RcppExport SEXP _STplm_logProfileCpp(SEXP theta0SEXP, SEXP CovStructureSEXP, SEXP loctimSEXP, SEXP DSSEXP, SEXP DTSEXP, SEXP XprimeSEXP, SEXP ZprimeSEXP, SEXP truncateSEXP, SEXP lambdaSEXP, SEXP ptermSEXP, SEXP DDnewSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type theta0(theta0SEXP);
    Rcpp::traits::input_parameter< const int >::type CovStructure(CovStructureSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type loctim(loctimSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type DS(DSSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type DT(DTSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type Xprime(XprimeSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type Zprime(ZprimeSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type truncate(truncateSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const int >::type pterm(ptermSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type DDnew(DDnewSEXP);
    rcpp_result_gen = Rcpp::wrap(logProfileCpp(theta0, CovStructure, loctim, DS, DT, Xprime, Zprime, truncate, lambda, pterm, DDnew));
    return rcpp_result_gen;
END_RCPP
}
// nonparlogprof
Eigen::MatrixXd nonparlogprof(const Eigen::VectorXd theta0, const int CovStructure, const Eigen::MatrixXd loctim, const Eigen::MatrixXd DS, const Eigen::MatrixXd DT, const Eigen::VectorXd Zprime, const Eigen::VectorXd truncate, const double lambda, const int pterm, const Eigen::VectorXd DDnew);
RcppExport SEXP _STplm_nonparlogprof(SEXP theta0SEXP, SEXP CovStructureSEXP, SEXP loctimSEXP, SEXP DSSEXP, SEXP DTSEXP, SEXP ZprimeSEXP, SEXP truncateSEXP, SEXP lambdaSEXP, SEXP ptermSEXP, SEXP DDnewSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type theta0(theta0SEXP);
    Rcpp::traits::input_parameter< const int >::type CovStructure(CovStructureSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type loctim(loctimSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type DS(DSSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type DT(DTSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type Zprime(ZprimeSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type truncate(truncateSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const int >::type pterm(ptermSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type DDnew(DDnewSEXP);
    rcpp_result_gen = Rcpp::wrap(nonparlogprof(theta0, CovStructure, loctim, DS, DT, Zprime, truncate, lambda, pterm, DDnew));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_STplm_covNonSep", (DL_FUNC) &_STplm_covNonSep, 5},
    {"_STplm_covNonSep2", (DL_FUNC) &_STplm_covNonSep2, 5},
    {"_STplm_covNonSep3", (DL_FUNC) &_STplm_covNonSep3, 6},
    {"_STplm_covSep", (DL_FUNC) &_STplm_covSep, 5},
    {"_STplm_loglikeCpp", (DL_FUNC) &_STplm_loglikeCpp, 10},
    {"_STplm_logProfileCpp", (DL_FUNC) &_STplm_logProfileCpp, 11},
    {"_STplm_nonparlogprof", (DL_FUNC) &_STplm_nonparlogprof, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_STplm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
