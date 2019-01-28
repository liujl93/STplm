#' Spatio-Temporal Partial Linear Model and Bandwidth Selection
#'
#' y(s,t) = X(s,t)'beta + f(t) + eps(s,t)
#' @param training a list containing space-time coordinates, design matrix X, response Z, distance in space Dloc, distance in time Dtim
#' @param CovStructure the covariance structure, 1 - a separable covariance matrix, 2 - a nonseparable covariance matrix
#' @param theta0 initial covariance paramters (not including sigma2)
#' @param lower,upper lower and upper limit of optimization for correlation parameters
#' @param lambda penalty coefficient
#' @param pterm number of parameters in penalty term
#' @param theta.true true correlation parameters
#' @param sigma2.true true sigma2
#' @return a list
#' \item{beta}{estimated regression coefficients}
#' \item{theta}{estimated covariance coefficients}
#' \item{J0betainv0}{...}
#' \item{J0theta}{...}
#' \item{Gaminv}{inverse of estimated covariance matrix * sigma2}
#' \item{likelihood}{negative log likelihood}
#' \item{convergence}{0-converge; 1-not converge; 2-theta known}
#' \item{CovStructure}{Covariance Structure}
#' @export

slr.fit = function(training,
                   CovStructure=2,theta0,lower=NULL,upper=NULL,
                   lambda=0,pterm=1,truncate.t=1,
                   theta.true=NULL,sigma2.true = NULL,DDnew=1) {
  loctim = as.matrix(training$loctim)
  P = nrow(loctim)
  X = training$X
  Z = training$Z
  Dloc = training$Dloc
  Dtim = training$Dtim
  
  if(is.null(theta.true)){
    fit = nlminb(theta0,logProfileCpp,lower=lower,upper=upper,
                 CovStructure=CovStructure,loctim=loctim,DS=Dloc,
                 DT=Dtim,Xprime=X,Zprime=Z,truncate=truncate.t,lambda=lambda,pterm=pterm,DDnew=DDnew)
    opt = fit$par; obj = fit$objective
    convergence = fit$convergence
  }else{
    opt = theta.true
    obj = logProfileCpp(theta.true, CovStructure=CovStructure,loctim=loctim,
                        DS=Dloc,DT=Dtim,Xprime=X,Zprime=Z,truncate=truncate.t,lambda=lambda,pterm=pterm,DDnew=DDnew)
    convergence = 2
  }

  ret = slr.est(training,CovStructure,opt,sigma2.true)
  ret$likelihood = obj
  ret$convergence = convergence
  ret$CovStructure = CovStructure
  return(ret)
  }

