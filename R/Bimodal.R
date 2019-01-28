#' @import locpol
#' Bimodal Kernel Function
#' Calculate the value of a bimodal kernel function
#' @param u a value
#' @return the value of bimodal kernel function at u
#' @export
K2 = function(u) {
  2/sqrt(pi)*u^2*exp(-u^2)
}
attr(K2,"dom") = c(-5,5)
attr(K2,"RK") = computeRK(K2)
attr(K2,"mu0K") = computeMu0(K2)
attr(K2,"mu2K") = computeMu(2,K2)
attr(K2,"K4") = computeK4(K2)

