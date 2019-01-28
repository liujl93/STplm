#' Calculate the Profile Likelihood
#'
#' Calculate the Profile Likelihood based on initial covariance parameters (not including sigma2)
#' @param theta0 initial covariance paramters (not including sigma2)
#' @param CovStructure the covariance structure, 1 - a separable covariance matrix, 2 - a nonseparable covariance matrix
#' @param Ds the distance matrix of location
#' @param Dt the distance matrix of time
#' @param P the dimension of data
#' @param Xprime transformed X, Xprime = (I-S)X
#' @param Zprime transformed Z, Zprime = (I-S)Z
#' @return the value of log profile likelihood.
#' @export
log.profile = function(theta0,CovStructure,loctim,Ds,Dt,Xprime,Zprime,Sigma=NULL) {
  P = nrow(Xprime)
  if (CovStructure == 1) {
    psi = (1-theta0[1])*exp(-Ds/theta0[2]-Dt/theta0[3])
    diag(psi) = 1 }
  else if (CovStructure == 2) {
    psi = (1-theta0[1])/(theta0[2]*(Dt)^2 + 1)*exp(- theta0[3]*Ds^2/(theta0[2]*(Dt)^2 + 1))
    diag(psi) = 1 }
  else if (CovStructure == 3) {
    psi <-  (1-theta0[1])/(theta0[2]^2*(Dt)^2+1)^(3/2)*exp(-theta0[3]*Ds)
    diag(psi) = 1 }
  else if (CovStructure == 4) {
    DD1 = theta0[4]*loctim[,3] + 1
    psi <-  DD1%*%t(DD1)/(theta0[2]^2*(Dt)^2+1)^(3/2)*exp(-theta0[3]*Ds)
    diag(psi) = diag(DD1%*% t(DD1))+theta0[1]  }
  else if (CovStructure == 5) {
    DD1 = theta0[4]*loctim[,3]+theta0[5]*loctim[,1]+theta0[6]*loctim[,2] + 1
    psi <-  DD1%*%t(DD1)/(theta0[2]^2*(Dt)^2+1)^(3/2)*exp(-theta0[3]*Ds)
    diag(psi) = diag(DD1%*% t(DD1))+theta0[1]   }
  else if (CovStructure == 6) {
    DD1 =1 + theta0[4]*loctim[,3]+theta0[5]*loctim[,3]^2+
      theta0[6]*loctim[,3]^3+theta0[7]*(ifelse((loctim[,3]-truncate.t)>0,(loctim[,3]-truncate.t),0))^3
    psi <-  DD1%*%t(DD1)/(theta0[2]^2*(Dt)^2+1)^(3/2)*exp(-theta0[3]*Ds)
    diag(psi) = diag(DD1%*% t(DD1)) +theta0[1]  }
  else stop("Please provide a valid covariance structure identification code (1-3).")
 U = base::chol(psi)
 sx = solve(t(U),Xprime)
 sz = solve(t(U),Zprime)
 beta = solve(t(sx)%*%sx,t(sx)%*%sz)
 ress = (sz - sx%*%beta)
 sigma2 = t(ress)%*%ress/P
 log.prof = P/2*log(sigma2) + sum(log(diag(U))) + P/2
 return(log.prof)
#  U <- chol(psi)
#  U.inv <- backsolve(U, diag(1, nrow = P))
#  psi.inv <- U.inv %*% t(U.inv)
#  # Compute the MLE for beta.
#  beta <- solve(t(Xprime) %*% psi.inv %*% Xprime) %*% t(Xprime) %*% psi.inv %*% Zprime
#  # Compute the MLE for sigma.
#  resid <- Zprime - Xprime %*% beta
#  sigma2 <- t(resid) %*% psi.inv %*% resid/P
#  # Evaluate -log-profile likelihood.
#  log.U.det <- sum(log(diag(U)))
#  # Log of the determinant of U.
#  prof <- log.U.det + (P * log(sigma2))/2 + P/2
#  return(prof)
}
