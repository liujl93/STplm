#' Fitting Spatio-Temporal Model with A Specified Bandwidth
#'
#' Fit a spatio-temporal model with a specified bandwidth
#' @param training a list containing space-time coordinates, design matrix X, response Z, distance in space Dloc, distance in time Dtim
#' @param theta0 initial covariance paramters (not including sigma2)
#' @param CovStructure the covariance structure, 1 - a separable covariance matrix, 2 - a nonseparable covariance matrix
#' @param h the prespecified bandwidth
#' @param KERN kernel function
#' @param lower,upper lower and upper limit of optimization for correlation parameters
#' @param lambda penalty coefficient
#' @param pterm number of parameters in penalty term
#' @param theta.true true correlation paramters, if provided, return theta
#' @param sigma2.true true sigma2, if theta.true is given, sigma2 must be given
#' @return beta.hat, theta.hat, f.hat, Sigma.hat, profile likelihood
#' @export
nonpar.spl = function(training,
                h,KERN,
                CovStructure,theta0,lower,upper,
                lambda=0,pterm=1,truncate.t=1,
                theta.true=NULL,sigma2.true=NULL,DDnew=1
                ) {
  loctim = as.matrix(training$loctim)
  P = nrow(loctim)
  X = training$X
  Z = training$Z
  Dloc = training$Dloc
  Dtim = training$Dtim
  
  # -------------------------------------- #    
  #      Smoother Matrix Calculation       #
  # -------------------------------------- #
  Smat = Smat1 = matrix(0,P,P)
    for(i in 1:P) {
      tmp = St(loctim[i,3],h,loctim[,3],KERN)
      Smat[i,] = tmp[1,]
      Smat1[i,] = tmp[2,]
    }
  
  # -------------------------------------- #    
  #            (I-S)X and (I-S)Z           #
  # -------------------------------------- #   
  Zprime = (diag(P)-Smat)%*%Z
  Xprime = (diag(P)-Smat)%*%X
  
  # -------------------------------------- #    
  #            Optimization                #
  # -------------------------------------- #
  if(is.null(theta.true)){
    # if theta.true is not provided, estimate using profile likelihood estimation
      fit = nlminb(theta0,nonparlogprof,lower=lower,upper=upper,
                   CovStructure=CovStructure,loctim=loctim,DS=Dloc,
                   DT=Dtim,Zprime=Zprime,truncate=truncate.t,lambda=lambda,pterm=pterm,DDnew=DDnew)
      opt = fit$par; obj = fit$objective
      convergence = fit$convergence
  }else{
    # if theta.true is provided, return theta and log-profile likelihood with theta
      opt = theta.true
      obj = nonparlogprof(theta.true, CovStructure=CovStructure,loctim=loctim,
                          DS=Dloc,DT=Dtim,Zprime=Zprime,truncate=truncate.t,lambda=lambda,pterm=pterm,DDnew=DDnew)
      convergence = 2 
  }
  # -------------------------------------- #    
  #                 \hat{theta}            #
  # -------------------------------------- #  
  results = nonpar.est(training,h,KERN,CovStructure,opt,sigma2.true)
  # -------------------------------------- #    
  #                 \hat{f,f1}             #
  # -------------------------------------- #   
  # fest = Smat%*%Z
  # fest1 = Smat1%*%Z/h
  # 
  # results$f.hat = fest
  # results$f1.hat = fest1
  return(results)
}
