#' Calculate the Smoothing Matrix
#'
#' Calculate the Smoothing Matrix
#' @param t a new time position
#' @param h bandwidth
#' @param FixT fixed time positions
#' @param KN Gaussian kernel or bimodal kernel
#' @return the smoothing matrix
#' @export
St = function(t,h,FixT,KN) {
  Ab = t(Dt(t,h,FixT))%*%Kt(t,h,FixT,KN)
  wt = solve(Ab%*%Dt(t,h,FixT),Ab)
  return(wt)
}
