#' Calculate Dt Matrix
#'
#' Calculate Dt matrix to be used in calculating smoothing matrix.
#' @param t a new time position
#' @param h bandwidth
#' @param FixT fixed time positions
#' @return a matrix
#' @export
Dt = function(t,h,FixT) cbind(rep(1,length(FixT)),(FixT-t)/h)
