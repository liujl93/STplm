#' @export
calc_ase <- function(time,ftrue,training,Type,object){
  if(Type %in% c(1,3)) KERN = "GK"
  if(Type %in% c(2,4)) KERN = "K2" 
  h <- object$bandwidth
  beta.hat <- object$beta.hat
  P0 = length(time)
  loctim = as.matrix(training$loctim)
  P = nrow(loctim)
  X = as.matrix(training$X)
  Z = as.numeric(training$Z)
  Dloc = training$Dloc
  Dtim = training$Dtim
  Smat = matrix(0,P0,P)
  for(i in 1:P0) {
    Smat[i,] = St(time[i],h,loctim[,3],KERN)[1,]
  }
  Zstar <- Z - X%*%beta.hat 
  fhat <- Smat%*%Zstar
  ase <- mean((fhat-ftrue)^2)
  return(ase)
}
