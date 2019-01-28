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
bw.cv = function(training,
                 h,KERN,
                 beta.fixed,
                 CovStructure,
                 theta.gcvce,Sigma.gcvce,
                 theta.gcvc=NULL,Sigma.gcvc=NULL
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

  theta.true = theta.gcvc[-1]
  theta.best = theta.gcvce[-1]

  if (CovStructure %in% 1:3) {
    DD0 = DD1 = rep(1,P)
  }else if (CovStructure == 4) {
    DD0 = theta.true[4]*loctim[,3]+1
    DD1 = theta.best[4]*loctim[,3]+1
  }else if (CovStructure == 5) {
    DD0 = theta.true[4]*loctim[,3]+theta.true[5]*loctim[,1]+theta.true[6]*loctim[,2]+1
    DD1 = theta.best[4]*loctim[,3]+theta.best[5]*loctim[,1]+theta.best[6]*loctim[,2]+1
  }else if (CovStructure == 6) {
    DD0 = theta.true[4]*loctim[,3]+theta.true[5]*(ifelse((loctim[,3]-truncate.t)>0,(loctim[,3]-truncate.t),0))+1
    DD1 = theta.best[4]*loctim[,3]+theta.best[5]*(ifelse((loctim[,3]-truncate.t)>0,(loctim[,3]-truncate.t),0))+1
  }else if (CovStructure == 7) {
    DD0 = theta.true[4]*loctim[,3]+theta.true[5]*loctim[,3]^2+theta.true[6]*(ifelse((loctim[,3]-truncate.t)>0,(loctim[,3]-truncate.t)^2,0))+1
    DD1 = theta.best[4]*loctim[,3]+theta.best[5]*loctim[,3]^2+theta.best[6]*(ifelse((loctim[,3]-truncate.t)>0,(loctim[,3]-truncate.t)^2,0))+1
  }else if (CovStructure == 8) {
    DD0 = theta.true[4]*loctim[,3]+theta.true[5]*loctim[,3]^2+theta.true[6]*loctim[,3]^3+theta.true[7]*(ifelse((loctim[,3]-truncate.t)>0,(loctim[,3]-truncate.t)^3,0))+1
    DD1 = theta.best[4]*loctim[,3]+theta.best[5]*loctim[,3]^2+theta.best[6]*loctim[,3]^3+theta.best[7]*(ifelse((loctim[,3]-truncate.t)>0,(loctim[,3]-truncate.t)^3,0))+1
  }else if (CovStructure == 9) {
    DD0 = theta.true[4]*(ifelse((loctim[,3]-truncate.t)>0,1,0))+1
    DD1 = theta.best[4]*(ifelse((loctim[,3]-truncate.t)>0,1,0))+1
  }  else stop("Invalid CovStructure (should be from 1 to 6).")

  if(is.null(Sigma.gcvc))  {GCVc = NA}
  else {
    div = theta.gcvc[1]*(mean(DD0^2)+theta.gcvc[2])
    GCVc  = mean((eprime/(1-sum(Smat*t(Sigma.gcvc)/div/P)))^2)
  }
  div = theta.gcvce[1]*(mean(DD1^2)+theta.gcvce[2])
  GCVce = mean((eprime/(1-sum(Smat*t(Sigma.gcvce)/div/P)))^2)
  
  result = list()
  result$CV = CV
  result$GCV = GCV
  result$GCVce = GCVce
  result$GCVc = GCVc

  return(result)
}
