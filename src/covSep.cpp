#include <RcppEigen.h>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <cmath>
#include <valarray>
#include <iostream>
#include "covSep.h"

using namespace Eigen;
using namespace Rcpp;

//' @export
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd covSep(const Eigen::MatrixXd DS, const Eigen::MatrixXd DT,
                       double c, double thetas, double thetat ){
    int P = DS.rows();
    Eigen::MatrixXd psi = (1-c)*(- DS.array()/thetas - DT.array()/thetat).exp().matrix();
    psi.diagonal() = Eigen::VectorXd(P).setOnes();
    return(psi);
}

/*** R
 # tests
 test <- covSep(A <- matrix(c(0,sqrt(2),sqrt(2),0),2), B <- matrix(c(0,sqrt(2),sqrt(2),0), ncol=2), 0.2, 0.3, 0.4)
 testR = (1-0.2)*exp(-A^2/0.3-(B*2)^2/0.4)
 diag(testR) = 1
 all.equal(test, testR)
 */
