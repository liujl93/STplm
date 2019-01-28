#' Spatio-Temporal Partial Linear Model and Bandwidth Selection
#'
#' y(s,t) = X(s,t)'beta + f(t) + eps(s,t)
#' @param training a list containing space-time coordinates, design matrix X, response Z, distance in space Dloc, distance in time Dtim
#' @param theta0 initial covariance paramters (not including sigma2)
#' @param CovStructure the covariance structure, 1 - a separable covariance matrix, 2 - a nonseparable covariance matrix
#' @param Type Scenario, 1 - GK + GK; 2 - K2 + K2; 3 - K2 + GK; 4 - GK + K2;
#' @param h0,h1 lower and upper limit of bandwidth
#' @param nlam number of bandwidth selection grid
#' @param lambda penalty coefficient
#' @param pterm number of parameters in penalty term
#' @param Sigma The true covariance matrix
#' @param lower,upper lower and upper limit of optimization for correlation parameters
#' @param theta.true true correlation parameters
#' @param sigma2.true true sigma2
#' @param beta.true true regression coefficients, used in step1 of bandwidth selection, oracle case
#' @return GCVce in detail only, CV samely.
#' \item{gcvce.bandwidth}{optimal bandwidth selected by GCVce criterion}
#' \item{gcvce.beta.est}{regression paramters estimator with bandwidth selected by GCVce criterion }
#' \item{gcvce.theta.est}{covariance paramters estimator with bandwidth selected by GCVce criterion }
#' \item{gcvce.f.est}{f(t) estimation at all t's with bandwidth selected by GCVce criterion}
#' \item{gcvce.likelihood}{log profile likelihood with bandwidth selected GCVce criterion}
#' @export

plm.cv = function(training,
                   h0=0.005,h1=0.1,nlam=30,Type=2,
                   CovStructure=2,theta0,lower=NULL,upper=NULL,
                   lambda=0,pterm=1,truncate.t=1,
                   theta.true=NULL,sigma2.true = NULL,
                   beta.true=NULL,DDnew=1) {
  loctim = as.matrix(training$loctim)
  P = nrow(loctim)
  X = as.matrix(training$X)
  Z = as.numeric(training$Z)
  Dloc = as.matrix(training$Dloc)
  Dtim = as.matrix(training$Dtim)
  
  if (Type == 1)       { KERN1 = "GK" ;KERN2 = "GK"
  }else if (Type == 2) { KERN1 = "K2" ;KERN2 = "K2"
  }else if (Type == 3) { KERN1 = "K2" ;KERN2 = "GK"
  }else if (Type == 4) { KERN1 = "GK" ;KERN2 = "K2"
  }else stop("Type Not Defined.")
  
  bw = exp(seq(log(h0),log(h1),length=nlam))
  cv=gcv=numeric(nlam)
  
  # ---------------------------------------------- #
  #               Bandwidth Selection              #
  # ---------------------------------------------- #
  
  # ----------- Step I (Use Kernel 1)  ----------- #
  
  # ---------------------------------------------- #
  #                 \tilde{beta}                   #
  # ---------------------------------------------- #  
  if(is.null(beta.true)){
    # Given h0, estimate \tilde{beta}
    # if theta is given
    result.ori = spl(training,
                     h0,KERN1,
                     CovStructure,theta0,lower,upper,
                     lambda,pterm,truncate.t,
                     theta.true=NULL,sigma2.true = NULL,
                     estimation=FALSE,DDnew=DDnew)
    beta.ori   = result.ori$beta.hat 
  }else{
    beta.ori = beta.true
  }
  
  # ----------- Step II (Use Kernel 2)  ---------- #
  
  # ---------------------------------------------- #
  #       y - x\tilde{beta} = f(t) + eps           #
  # ---------------------------------------------- #   
  for(i in 1:nlam) {
    training0 = training
    training0$Z = Z - X%*%beta.ori
    bw.fit = bw.cvonly(training,
                   bw[i],KERN2,
                   beta.ori)
    cv[i] = bw.fit$CV
    gcv[i] = bw.fit$GCV
  }
  
  # ----------- Step III (Use Kernel 2)  --------- #
  
  
  ret = list()
  ret$bw = cbind(cv,gcv)

  hopt <- bw[apply(ret$bw,2,which.min)]
  hoptuni <- unique(hopt)
  
  fit = list()
  for(iter in 1:length(hoptuni)){
    fit[[iter]]=spl(training,
                    hoptuni[iter],KERN2,
                    CovStructure,theta0,lower,upper,
                    lambda,pterm,truncate.t,
                    theta.true,sigma2.true,
                    estimation=TRUE,DDnew=DDnew) 
    fit[[iter]]$bandwidth = hoptuni[iter]
    fit[[iter]]$CovStructure = CovStructure
  }
  # ------------------- CV ------------------- #
  tmp = which(hoptuni == hopt[1])
  fit1 = fit[[tmp]]
  # ------------------- GCV ------------------ #
  tmp = which(hoptuni == hopt[2])
  fit2 = fit[[tmp]]

  ret$cv = fit1
  ret$gcv = fit2

  return(ret)
}
