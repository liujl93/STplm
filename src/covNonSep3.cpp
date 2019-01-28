#include <RcppEigen.h>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <cmath>
#include <valarray>
#include <iostream>
#include "covNonSep2.h"

using namespace Eigen;
using namespace Rcpp;

//' @export
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd covNonSep3(const Eigen::MatrixXd DS, const Eigen::MatrixXd DT,
                          double c, double thetat, double thetas, double theta3){
    int P = DS.rows();
    Eigen::MatrixXd psi = ( (1-c)*theta3 / (((thetat*DT.array()).pow(2) + 1).sqrt()) / ((thetat*DT.array()).pow(2) + theta3) * (- thetas * (((thetat*DT.array()).pow(2) + 1)/((thetat* DT.array()).pow(2) + theta3)).sqrt() * DS.array()).exp() ).matrix();
    psi.diagonal() = Eigen::VectorXd(P).setOnes();
    return(psi);
}

/*** R
 # tests
 c = 0.2
 thetas = 0.3
 thetat = 0.4
 theta3 = 4
 test <- covNonSep3(A <- matrix(c(0,sqrt(2),sqrt(2),0),2), B <- matrix(c(0,sqrt(2),sqrt(2),0), ncol=2), c,  thetat, thetas, theta3)
 testR = (1-c)*theta3/(thetat^2*B^2+1)^(1/2)/(thetat^2*B^2+theta3)*exp(-thetas*sqrt((thetat^2*B^2+1)/(thetat^2*B^2+theta3))*A)
 diag(testR) = 1
 all.equal(test, testR)
 */
