#' Estimate beta, sigma2, f(t) and covariance matrix
#'
#' Estimate beta, sigma2, f(t) and covariance matrix after obtaining the profile likelihood estimator of theta (sigma2 excluded).
#' @param theta.best profile likelihood estimator of theta (not including sigma2)
#' @param CovStructure specify the covariance structure, 1 - a separable covariance matrix, 2 - a nonseparable covariance matrix
#' @param Ds the distance matrix of location
#' @param Dt the distance matrix of time
#' @param P the dimension of data
#' @param X the design matrix of linear part of the model
#' @param Z the dependent variable
#' @keywords profile likelihood
#' @return
#' \item{beta}{Linear regression paramters estimation}
#' \item{theta}{Covariance paramters estimation}
#' \item{fest}{f(t) estimation}
#' \item{Sigma}{Estimated covariance matrix}
#' \item{Gaminv}
#' @import compiler
#' @export

nonpar.est = function(training,
                      h,KERN,
                      CovStructure,theta.best,
                      sigma2.true=NULL) {
  # ---------------------------------------------- #    
  #                 Y = f(t) + eps                 #
  #     Method: Profile Likelihood Estimation      #
  # ---------------------------------------------- #   
  loctim = as.matrix(training$loctim)
  P = nrow(loctim)
  Ds = training$Dloc
  Dt = training$Dtim
  Z = training$Z
  
  # ---------------------------------------------- #    
  #  Calculate correlation matrix (Gamma/sigma^2)  #
  # ---------------------------------------------- #   
  if (CovStructure == 1) {
    DD1 = 1
    psi = (1-theta.best[1])*exp(-Ds/theta.best[2]-Dt/theta.best[3])
    diag(psi) = 1 
  }else if (CovStructure == 2) {
    DD1 = 1
    psi <-  theta.best[4]/(theta.best[2]^2*Dt^2+1)^(1/2)/(theta.best[2]^2*Dt^2+theta.best[4])*exp(-theta.best[3]*Ds*((theta.best[2]^2*Dt^2+1)/(theta.best[2]^2*Dt^2+theta.best[4]))^(1/2))
    diag(psi) = 1 + theta.best[1]
  }else if (CovStructure == 3) {
    DD1 = 1
    psi <-  (1-theta.best[1])/(theta.best[2]^2*Dt^2+1)^(3/2)*exp(-theta.best[3]*Ds)
    diag(psi) = 1 
  }else if (CovStructure == 4) {
    DD1 = theta.best[4]*loctim[,3]+1
    psi <-  DD1%*%t(DD1)/(theta.best[2]^2*Dt^2+1)^(3/2)*exp(-theta.best[3]*Ds)
    diag(psi) = diag(DD1%*% t(DD1))+theta.best[1]
  }else if (CovStructure == 5) {
    DD1 = theta.best[4]*loctim[,3]+theta.best[5]*loctim[,1]+theta.best[6]*loctim[,2]+1
    psi <-  DD1%*%t(DD1)/(theta.best[2]^2*Dt^2+1)^(3/2)*exp(-theta.best[3]*Ds)
    diag(psi) = diag(DD1%*% t(DD1))+theta.best[1]
  }else if (CovStructure == 6) {
    DD1 = theta.best[4]*loctim[,3]+theta.best[5]*(ifelse((loctim[,3]-truncate.t)>0,(loctim[,3]-truncate.t),0))+1
    psi <-  DD1%*%t(DD1)/(theta.best[2]^2*Dt^2+1)^(3/2)*exp(-theta.best[3]*Ds)
    diag(psi) = diag(DD1%*% t(DD1))+theta.best[1]
  }else if (CovStructure == 7) {
    DD1 = theta.best[4]*loctim[,3]+theta.best[5]*loctim[,3]^2+theta.best[6]*(ifelse((loctim[,3]-truncate.t)>0,(loctim[,3]-truncate.t)^2,0))+1
    psi <-  DD1%*%t(DD1)/(theta.best[2]^2*Dt^2+1)^(3/2)*exp(-theta.best[3]*Ds)
    diag(psi) = diag(DD1%*% t(DD1))+theta.best[1]
  }else if (CovStructure == 8) {
    DD1 = theta.best[4]*loctim[,3]+theta.best[5]*loctim[,3]^2+theta.best[6]*loctim[,3]^3+theta.best[7]*(ifelse((loctim[,3]-truncate.t)>0,(loctim[,3]-truncate.t)^3,0))+1
    psi <-  DD1%*%t(DD1)/(theta.best[2]^2*Dt^2+1)^(3/2)*exp(-theta.best[3]*Ds)
    diag(psi) = diag(DD1%*% t(DD1))+theta.best[1]
  }  else stop("Invalid CovStructure (should be from 1 to 6).")

  # ----------------------------------------------------- #    
  #    Inverse correlation matrix (Gamma/sigma^2)^(-1)    #
  # ----------------------------------------------------- #      
  U = base::chol(psi)
  U.inv = backsolve(U,diag(1,nrow=P))
  Gaminv = U.inv%*%t(U.inv)
  
  # -------------------------------------- #    
  #                 \hat{sigma2}           #
  # -------------------------------------- #  
  if(is.null(sigma2.true)){
      sigma2 = as.numeric(t(Z)%*%Gaminv%*%(Z)/P)
  }
  else{
    sigma2 = sigma2.true
  }
  
  Sigma = sigma2*psi

  result = list()
  result$theta.hat = c(sigma2,theta.best)
  result$Sigma.hat = Sigma
  #result$Gaminv = Gaminv
  return(result)
}

