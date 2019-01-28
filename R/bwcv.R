#' Bandwidth Selection
#'
#' Bandwidth Selection Using Cross-Valication and Generalized Cross-Validation Correlated Criterion
#' @param h the prespecified bandwidth
#' @param KERN kernel function
#' @param loctim (s,t) the spatiotemporal index
#' @param P the dimension of data
#' @param X X
#' @param Z y
#' @param beta.fixed,Sigma.fixed,sigma2.fixed prespecified parameters and covariance matrix
#' @return CV, GCVce criterion
#' @import compiler
#' @export
bw.cvonly = function(training,
                 h,KERN,
                 beta.fixed
) {
  loctim = as.matrix(training$loctim)
  P = nrow(loctim)
  X = training$X
  Z = training$Z
  
  # -------------------------------------- #
  #      Smoother Matrix Calculation       #
  # -------------------------------------- #
  Smat = matrix(0,P,P)
  for(i in 1:P) {
    Smat[i,] = St(loctim[i,3],h,loctim[,3],KERN)[1,]
  }
  
  # -------------------------------------- #
  #            (I-S)X and (I-S)Z           #
  # -------------------------------------- #
  Zprime = (diag(P)-Smat)%*%Z
  Xprime = (diag(P)-Smat)%*%X
  
  eprime = Zprime - Xprime%*%beta.fixed
  
  CV    = mean((eprime/(1-diag(Smat)))^2)
  GCV   = mean((eprime/(1-1/P*sum(diag(Smat))))^2)
  

  result = list()
  result$CV = CV
  result$GCV = GCV

  return(result)
}
