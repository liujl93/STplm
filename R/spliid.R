#' @import compiler
#' @export
spl.iid <- function(training, h, KERN) {
  loctim = as.matrix(training$loctim)
  P = nrow(loctim)
  X = training$X
  Z = training$Z
  
  Smat = Smat1 = matrix(0,P,P)
  for(i in 1:P) {
    tmp = St(loctim[i,3],h,loctim[,3],KERN)
    Smat[i,] = tmp[1,]
    Smat1[i,] = tmp[2,]
  }

  Zprime = (diag(P) - Smat) %*% Z
  Xprime = (diag(P) - Smat) %*% X
  
  beta.est = solve(t(Xprime) %*% Xprime, t(Xprime) %*% Zprime)
  sigma2 = as.numeric(t(Zprime - Xprime %*% beta.est) %*% (Zprime - Xprime %*% beta.est)/P)
  fest = Smat %*% (Z - X %*% beta.est)
  fest1 = Smat1 %*% (Z - X %*% beta.est)/h
  
  J0betainv0 = solve(t(Xprime)%*%Xprime)*sigma2 # X(I-S)Gaminv(I-S)X

  
  J0f = numeric(P)
  psi = NULL
  
  J0f_qt = J0f_noqt = J0f_plugin = numeric(P)
  for (i in 1:P) {
    J0f_qt[i] = kkt(loctim[i,3],h,loctim[,3],KERN,psi,sigma2,method="qt")
    J0f_noqt[i] = kkt(loctim[i,3],h,loctim[,3],KERN,psi,sigma2,method="noqt")
    J0f_plugin[i] = kkt(loctim[i,3],h,loctim[,3],KERN,psi,sigma2,method="plugin")
  }
  
  ret = list()
  ret$beta.hat = beta.est
  ret$theta.hat = sigma2
  ret$f.hat = fest
  ret$f.hat1 = fest1
  ret$J0betainv0 = J0betainv0
  ret$J0f_qt = J0f_qt
  ret$J0f_noqt = J0f_noqt
  ret$J0f_plugin = J0f_plugin
  return(ret)
}


