#include <RcppEigen.h>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <cmath>
#include <valarray>
#include <iostream>
#include "covNonSep.h"

using namespace Eigen;
using namespace Rcpp;

//' @export
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd covNonSep(const Eigen::MatrixXd DS, const Eigen::MatrixXd DT,
                       double c, double thetas, double thetat ){
    int P = DS.rows();
    Eigen::MatrixXd psi = (( (1-c) / (thetat*(DT.array()).pow(2) + 1) ).cwiseProduct( (- thetas*DS.array().pow(2) / (thetat*(DT.array()).pow(2) + 1 ) ).exp())).matrix();
    psi.diagonal() = Eigen::VectorXd(P).setOnes();
    return(psi);
}

/*** R
 # tests
 theta0 = c(0.2,0.3,0.4)
 P = 6
 loctim = cbind(sx = c(1:6), sy = c(1:6), t = c(1:6))
 DS = as.matrix(dist(loctim[,-3]))
 DT = as.matrix(dist(loctim[,3]),p=1)
 test <- covNonSep(DS, DT, theta0[1],theta0[2],theta0[3])
 testR = (1-theta0[1])/(theta0[3]*(P*DT)^2+1)*exp(-theta0[2]*DS^2/(theta0[3]*(P*DT)^2+1))
 diag(testR) = 1
 sum(test != testR)
 */
