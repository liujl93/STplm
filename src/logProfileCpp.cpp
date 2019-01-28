 #include <RcppEigen.h>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <cmath>        // std::exp(double)
#include <valarray>     // std::valarray, std::exp(valarray)
#include <iostream>
#include <math.h>
#include <Rmath.h>
#include "covSep.h"
#include "covNonSep.h"
#include "covNonSep2.h"
#include "covNonSep3.h"

using namespace Eigen;
using namespace Rcpp;

//' @export
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd logProfileCpp(const Eigen::VectorXd theta0,
                     const int CovStructure,
                     const Eigen::MatrixXd loctim,
                     const Eigen::MatrixXd DS,
                     const Eigen::MatrixXd DT,
                     const Eigen::MatrixXd Xprime,
                     const Eigen::VectorXd Zprime,
                     const Eigen::VectorXd truncate,
                     const double lambda,
                     const int pterm,
                     const Eigen::VectorXd DDnew
                              ) {

    int P = Zprime.size();
    double Pd = Zprime.size();
    MatrixXd psi(MatrixXd(P,P).setZero());
    if(CovStructure == 1){
        psi = covSep(DS,DT,theta0(0),theta0(1),theta0(2));
    } else if(CovStructure == 2){
        psi = covNonSep3(DS,DT,0,theta0(1),theta0(2),theta0(3));
        psi.diagonal() = (1+theta0(0))*VectorXd(P).setOnes();
    } else if(CovStructure == 3) {
        psi = covNonSep2(DS,DT,0,theta0(1),theta0(2));
        psi.diagonal() = (1+theta0(0))*VectorXd(P).setOnes();
    } else if(CovStructure == 4) {
        Eigen::VectorXd DD1 = theta0(3)*loctim.col(2).array() + 1;
        psi = (DD1 * DD1.adjoint()).cwiseProduct(covNonSep2(DS,DT,0,theta0(1),theta0(2)));
        psi.diagonal() = (DD1 * DD1.adjoint()).diagonal() + theta0(0)*VectorXd(P).setOnes();
    } else if(CovStructure == 5) {
        Eigen::VectorXd DD1 = theta0(3)*loctim.col(2).array() + theta0(4)*loctim.col(0).array() + theta0(5)*loctim.col(1).array() + 1;
        psi = (DD1 * DD1.adjoint()).cwiseProduct(covNonSep2(DS,DT,0,theta0(1),theta0(2)));
        psi.diagonal() = (DD1 * DD1.adjoint()).diagonal()+ theta0(0)*VectorXd(P).setOnes();
    } else if(CovStructure == 6) {
        Eigen::VectorXd DD1 = 1 + theta0(3)*loctim.col(2).array() + theta0(4)*(loctim.col(2).array() > truncate(0)).select((loctim.col(2).array() - truncate(0)), 0);
        psi = (DD1 * DD1.adjoint()).cwiseProduct(covNonSep2(DS,DT,0,theta0(1),theta0(2)));
        psi.diagonal() = (DD1 * DD1.adjoint()).diagonal()+ theta0(0)*VectorXd(P).setOnes();
    } else if(CovStructure == 7) {
        Eigen::VectorXd DD1 = 1 + theta0(3)*loctim.col(2).array() + theta0(4)*loctim.col(2).array().pow(2) + theta0(5)*(loctim.col(2).array() > truncate(0)).select((loctim.col(2).array() - truncate(0)).pow(2), 0);
        psi = (DD1 * DD1.adjoint()).cwiseProduct(covNonSep2(DS,DT,0,theta0(1),theta0(2)));
        psi.diagonal() = (DD1 * DD1.adjoint()).diagonal()+ theta0(0)*VectorXd(P).setOnes();
    } else if(CovStructure == 8) {
        Eigen::VectorXd DD1 = 1 + theta0(3)*loctim.col(2).array() + theta0(4)*loctim.col(2).array().pow(2) + theta0(5)*loctim.col(2).array().pow(3) + theta0(6)*(loctim.col(2).array() > truncate(0)).select((loctim.col(2).array() - truncate(0)).pow(3), 0);
        psi = (DD1 * DD1.adjoint()).cwiseProduct(covNonSep2(DS,DT,0,theta0(1),theta0(2)));
        psi.diagonal() = (DD1 * DD1.adjoint()).diagonal()+ theta0(0)*VectorXd(P).setOnes();
    } else if(CovStructure == 9) {
        Eigen::VectorXd DD1 = 1 + theta0(3)*(loctim.col(2).array() > truncate(0)).select(VectorXd(P).setOnes().array(),0);
        psi = (DD1 * DD1.adjoint()).cwiseProduct(covNonSep2(DS,DT,0,theta0(1),theta0(2)));
        psi.diagonal() = (DD1 * DD1.adjoint()).diagonal()+ theta0(0)*VectorXd(P).setOnes();
    } else if(CovStructure == 10) {
        Eigen::VectorXd DD1 = 1 + theta0(3)*(loctim.col(2).array() > truncate(0)).select((loctim.col(2).array() - truncate(0)),0) + theta0(4)*(loctim.col(2).array() > truncate(1)).select((loctim.col(2).array() - truncate(1)),0);
        psi = (DD1 * DD1.adjoint()).cwiseProduct(covNonSep2(DS,DT,0,theta0(1),theta0(2)));
        psi.diagonal() = (DD1 * DD1.adjoint()).diagonal()+ theta0(0)*VectorXd(P).setOnes();
    } else if(CovStructure == 11) {
        psi = (DDnew * DDnew.adjoint()).cwiseProduct(covNonSep2(DS,DT,0,theta0(1),theta0(2)));
        psi.diagonal() = (DDnew * DDnew.adjoint()).diagonal()+ theta0(0)*VectorXd(P).setOnes();
    } else if(CovStructure == 12) {
        Eigen::VectorXd DD1 = 1 + theta0(3)*loctim.col(2).array() + theta0(4)*(loctim.col(2).array() > truncate(0)).select((loctim.col(2).array() - truncate(0)),0) + theta0(5)*(loctim.col(2).array() > truncate(1)).select((loctim.col(2).array() - truncate(1)),0);
        psi = (DD1 * DD1.adjoint()).cwiseProduct(covNonSep2(DS,DT,0,theta0(1),theta0(2)));
        psi.diagonal() = (DD1 * DD1.adjoint()).diagonal()+ theta0(0)*VectorXd(P).setOnes();
    }
    Eigen::MatrixXd U = psi.llt().matrixL().adjoint();
    Eigen::MatrixXd SX = psi.llt().solve(Xprime);
    Eigen::VectorXd beta = ((Xprime.adjoint())*SX).llt().solve((SX.adjoint())*Zprime);
    Eigen::VectorXd res = Zprime - Xprime*beta;
    double sigma2 = (res.adjoint())*(psi.llt().solve(res));
    sigma2 = sigma2/Pd;
    Eigen::MatrixXd logprof(Eigen::MatrixXd(1,1).setZero());
    logprof(0,0) = (Pd*log(sigma2)/2 + U.diagonal().array().log().sum() + Pd/2 + lambda*theta0.tail(pterm).norm());
    return(logprof);
}

/*** R
 # Test:
 library(STplm)
 theta0 = c(0.2,0.3,0.4)
 P = 6
 CovStructure = 1
 loctim = cbind(sx = c(1:6), sy = c(1:6), t = c(1:6))
 DS = as.matrix(dist(loctim[,-3]))
 DT = as.matrix(dist(loctim[,3]),p=1)
 Zprime = rnorm(6)
 Xprime = matrix(rnorm(6*3),6)
 system.time(test <- logProfileCpp(theta0, CovStructure,loctim,DS,DT,Xprime,Zprime,1,0,1))
 system.time(testR <- log.profile(theta0, CovStructure,loctim,DS,DT,Xprime,Zprime))
 test == testR
 */
