#' Calculate Kt Matrix in PLM
#'
#' Calculate Kt matrix to be used in calculating the smoothing matrix
#' @param t a new time position
#' @param h bandwidth
#' @param FixT fixed time positions
#' @param KN Gaussian kernel or bimodal kernel
#' @return a matrix
#' @export
Kt = function(t,h,FixT,KN) {
  switch(KN,
         GK = diag(dnorm((FixT-t)/h)/h),
         K2 = diag(K2((FixT-t)/h))/h)
}
