#' @import KernSmooth
#' @export
predict.zm <- function(object,training,test){
  # Fitted Model: object
  # iid case
  loctim = as.matrix(training$loctim)
  P = nrow(loctim)
  Z = as.numeric(training$Z)
  # Test data
  sp0 = as.matrix(test$loctim[,-ncol(test$loctim)])
  time0 = as.numeric(test$loctim[, ncol(test$loctim)])
  P0 = nrow(sp0)
  
  Dloc0 = proxy::dist(sp0, loctim[, -ncol(loctim)]) * l
  Dtim0 = proxy::dist(time0, loctim[, ncol(loctim)]) * constt
  

  if(is.null(object$CovStructure)){
    pred = 0
  }  else{
    CovStructure = object$CovStructure
    theta.hat = object$theta
    Gaminv.hat = object$Gaminv
    
    if (CovStructure == 1) {
      DD1 = 1
      Gam0 = (1 - theta.hat[2]) * exp(-Dloc0/theta.hat[3] - Dtim0/theta.hat[4])
    }
    if (CovStructure == 2) {
      DD1 = 1
      Gam0 = theta.hat[5]/(theta.hat[3]^2*Dtim0^2+1)^(1/2)/(theta.hat[3]^2*Dtim0^2+theta.hat[5])*exp(-theta.hat[4]*Dloc0*((theta.hat[3]^2*Dtim0^2+1)/(theta.hat[3]^2*Dtim0^2+theta.hat[5]))^(1/2))
    }
    if (CovStructure == 3) {
      DD1 = 1
      Gam0 <- (1-theta.hat[2])/(theta.hat[3]^2*(Dtim0)^2+1)^(3/2)*exp(-theta.hat[4]*Dloc0)
    }
    if (CovStructure == 4) {
      DD1 = theta.hat[5] * time0 + 1
      DD2 = theta.hat[5] * loctim[, 3] + 1
      Gam0 <- DD1 %*% t(DD2)/(theta.hat[3]^2*(Dtim0)^2 + 1)^(3/2) * exp(-theta.hat[4] * Dloc0)
    }
    if (CovStructure == 5) {
      DD1 = theta.hat[5]*time0 + theta.hat[6]*sp0[,1] + theta.hat[7]*sp0[,2] + 1
      DD2 = theta.hat[5]*loctim[,3] + theta.hat[6]*loctim[,1] + theta.hat[7]*loctim[,2] + 1
      Gam0 <- DD1 %*% t(DD2)/(theta.hat[3]^2*(Dtim0)^2 + 1)^(3/2) * exp(-theta.hat[4] * Dloc0)
    }
    if (CovStructure == 6) {
      DD1 =1 + theta.hat[7]*time0^3 + theta.hat[8]*(ifelse((time0-truncate.t)>0,(time0-truncate.t),0))^3
      DD2 =1 + theta.hat[7]*(loctim[,3])^3 + theta.hat[8]*(ifelse((loctim[,3]-truncate.t)>0,(loctim[,3]-truncate.t),0))^3
      Gam0 <- DD1 %*% t(DD2)/(theta.hat[3]^2*(Dtim0)^2 + 1)^(3/2) * exp(-theta.hat[4] * Dloc0)
    }
    
    pred = Gam0 %*% Gaminv.hat %*% Z
  }
  
  return(pred)
}
