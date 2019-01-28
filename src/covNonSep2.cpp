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
Eigen::MatrixXd covNonSep2(const Eigen::MatrixXd DS, const Eigen::MatrixXd DT,
                          double c, double thetat, double thetas ){
    int P = DS.rows();
    Eigen::MatrixXd psi = ((1-c)/(((thetat*DT.array()).pow(2) + 1).pow(1.5)) * (- thetas*DS.array()).exp()).matrix();
    psi.diagonal() = Eigen::VectorXd(P).setOnes();
    return(psi);
}

/*** R
 # tests
 c = 0.2
 thetas = 0.3
 thetat = 0.4
 test <- covNonSep2(A <- matrix(c(0,sqrt(2),sqrt(2),0),2), B <- matrix(c(0,sqrt(2),sqrt(2),0), ncol=2), 0.2, 0.3, 0.4)
 testR = (1-c)/((thetat*2*B)^2+1)^(3/2)*exp(-thetas*A)
 diag(testR) = 1
 all.equal(test, testR)
 */
