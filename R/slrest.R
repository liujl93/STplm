#' Estimate beta, sigma2, f(t) and covariance matrix
#'
#' Estimate beta, sigma2, f(t) and covariance matrix after obtaining the profile likelihood estimator of theta (sigma2 excluded).
#' @param training a list containing space-time coordinates, design matrix X, response Z, distance in space Dloc, distance in time Dtim
#' @param theta.best profile likelihood estimator of theta (not including sigma2)
#' @param CovStructure specify the covariance structure, 1 - a separable covariance matrix, 2 - a nonseparable covariance matrix
#' @param Sigma true covariance matrix 
#' @keywords profile likelihood
#' @return
#' \item{beta}{Linear regression paramters estimation}
#' \item{theta}{Covariance paramters estimation}
#' \item{fest}{f(t) estimation}
#' \item{Sigma}{Estimated covariance matrix}
#' \item{Gaminv}
#' @import compiler
#' @export

slr.est = function(training,
               CovStructure,theta.best,
               sigma2.true=NULL,all=FALSE) {
  loctim = as.matrix(training$loctim)
  P = nrow(loctim)
  X = training$X
  Z = training$Z
  Ds = training$Dloc
  Dt = training$Dtim
  
  if (CovStructure == 1) {
      DD1 = 1
      psi = (1-theta.best[1])*exp(-Ds/theta.best[2]-Dt/theta.best[3])
      diag(psi) = 1
  }
  if (CovStructure == 2) {
      DD1 = 1
      psi <-  theta.best[4]/(theta.best[2]^2*Dt^2+1)^(1/2)/(theta.best[2]^2*Dt^2+theta.best[4])*exp(-theta.best[3]*Ds*((theta.best[2]^2*Dt^2+1)/(theta.best[2]^2*Dt^2+theta.best[4]))^(1/2))
      diag(psi) = 1 + theta.best[1]
  }
  if (CovStructure == 3) {
      DD1 = 1
      psi <-  1/(theta.best[2]^2*Dt^2+1)^(3/2)*exp(-theta.best[3]*Ds)
      diag(psi) = 1 + theta.best[1]
  }
  if (CovStructure == 4) {
      DD1 = theta.best[4]*loctim[,3]+1
      psi <-  DD1%*%t(DD1)/(theta.best[2]^2*Dt^2+1)^(3/2)*exp(-theta.best[3]*Ds)
      diag(psi) = diag(DD1%*% t(DD1))+theta.best[1]
  }
  if (CovStructure == 5) {
      DD1 = theta.best[4]*loctim[,3]+theta.best[5]*loctim[,1]+theta.best[6]*loctim[,2]+1
      psi <-  DD1%*%t(DD1)/(theta.best[2]^2*Dt^2+1)^(3/2)*exp(-theta.best[3]*Ds)
      diag(psi) = diag(DD1%*% t(DD1))+theta.best[1]
  }
  if (CovStructure == 6) {
      DD1 = theta.best[4]*loctim[,3]+theta.best[5]*(ifelse((loctim[,3]-truncate.t)>0,(loctim[,3]-truncate.t),0))+1
      psi <-  DD1%*%t(DD1)/(theta.best[2]^2*Dt^2+1)^(3/2)*exp(-theta.best[3]*Ds)
      diag(psi) = diag(DD1%*% t(DD1))+theta.best[1]
  }
  if (CovStructure == 7) {
      DD1 = theta.best[4]*loctim[,3]+theta.best[5]*loctim[,3]^2+theta.best[6]*(ifelse((loctim[,3]-truncate.t)>0,(loctim[,3]-truncate.t)^2,0))+1
      psi <-  DD1%*%t(DD1)/(theta.best[2]^2*Dt^2+1)^(3/2)*exp(-theta.best[3]*Ds)
      diag(psi) = diag(DD1%*% t(DD1))+theta.best[1]
  }
  if (CovStructure == 8) {
      DD1 = theta.best[4]*loctim[,3]+theta.best[5]*loctim[,3]^2+theta.best[6]*loctim[,3]^3+theta.best[7]*(ifelse((loctim[,3]-truncate.t)>0,(loctim[,3]-truncate.t)^3,0))+1
      psi <-  DD1%*%t(DD1)/(theta.best[2]^2*Dt^2+1)^(3/2)*exp(-theta.best[3]*Ds)
      diag(psi) = diag(DD1%*% t(DD1))+theta.best[1]
  }
  if (CovStructure == 10) {
      DD1 = theta.best[4]*(ifelse((loctim[,3]-truncate.t[1])>0,(loctim[,3]-truncate.t[1]),0))+theta.best[5]*(ifelse((loctim[,3]-truncate.t[2])>0,(loctim[,3]-truncate.t[2]),0))+1
      psi <-  DD1%*%t(DD1)/(theta.best[2]^2*Dt^2+1)^(3/2)*exp(-theta.best[3]*Ds)
      diag(psi) = diag(DD1%*% t(DD1))+theta.best[1]
  }
  if (CovStructure == 11) {
      DD1 = DDnew
      psi <-  DD1%*%t(DD1)/(theta.best[2]^2*Dt^2+1)^(3/2)*exp(-theta.best[3]*Ds)
      diag(psi) = diag(DD1%*% t(DD1))+theta.best[1]
  }
  if (CovStructure == 12) {
      DD1 = theta.best[4]*loctim[,3]+theta.best[5]*(ifelse((loctim[,3]-truncate.t[1])>0,(loctim[,3]-truncate.t[1]),0))+theta.best[6]*(ifelse((loctim[,3]-truncate.t[2])>0,(loctim[,3]-truncate.t[2]),0))+1
      psi <-  DD1%*%t(DD1)/(theta.best[2]^2*Dt^2+1)^(3/2)*exp(-theta.best[3]*Ds)
      diag(psi) = diag(DD1%*% t(DD1))+theta.best[1]
  } 
  U = base::chol(psi)
  U.inv = backsolve(U,diag(1,nrow=P))
  
  Gaminv = U.inv%*%t(U.inv)

  beta.est = solve(t(X)%*%Gaminv%*%X,t(X)%*%Gaminv%*%Z)

  if(is.null(sigma2.true)){
      sigma2 = as.numeric(t(Z-X%*%beta.est)%*%Gaminv%*%(Z-X%*%beta.est)/P)
  }else{
      sigma2 = sigma2.true
  }
  
  J0betainv0 = solve(t(X)%*%Gaminv%*%X)*sigma2 
  
  # ------------------ I(theta) ------------------ # 
  if(CovStructure == 1){
    psic = -psi/(1-theta.best[1])
    diag(psic) = 0
    
    psia = psi*Ds/theta.best[2]^2
    psia[which(Ds==0,arr.ind = TRUE)] = 0
    
    psib = psi*Dt/theta.best[3]^2
    psib[which(Dt==0,arr.ind = TRUE)] = 0
    
    J0theta = matrix(0,4,4)
    g1 = (Gaminv%*%psic)
    g2 = (Gaminv%*%psia)
    g3 = (Gaminv%*%psib)
    
    J0theta[1,1] = sum(g1*t(g1))/2
    J0theta[2,2] = sum(g2*t(g2))/2
    J0theta[3,3] = sum(g3*t(g3))/2
    J0theta[4,4] = P/2/sigma2^2
    
    J0theta[1,2] = J0theta[2,1] = sum(g1*t(g2))/2
    J0theta[1,3] = J0theta[3,1] = sum(g1*t(g3))/2
    J0theta[2,3] = J0theta[3,2] = sum(g2*t(g3))/2
    J0theta[1,4] = J0theta[4,1] = sum(diag(1/sigma2,P)*t(g1))/2
    J0theta[2,4] = J0theta[4,2] = sum(diag(1/sigma2,P)*t(g2))/2
    J0theta[3,4] = J0theta[4,3] = sum(diag(1/sigma2,P)*t(g3))/2
    #J0thetainv = solve(J0theta)
  }
  if(CovStructure == 2){
    psic <- diag(P)
    psi1 <- -psi*theta.best[2]*Dt^2/(theta.best[2]^2 *(Dt)^2 + 1)/(theta.best[2]^2 *(Dt)^2 + theta.best[4])*(3*theta.best[2]^2*(Dt)^2+theta.best[4]+2+theta.best[3]*Ds*(theta.best[4]-1)*sqrt(theta.best[2]^2 *(Dt)^2 + 1)/sqrt(theta.best[2]^2 *(Dt)^2 + theta.best[4]))
    psi1[which(Dt==0,arr.ind = TRUE)] <- 0
    psi2 <- -psi*sqrt(theta.best[2]^2 *(Dt)^2 + 1)/sqrt(theta.best[2]^2 *(Dt)^2 + theta.best[4])*Ds
    psi2[which(Ds==0,arr.ind = TRUE)] <- 0
    psi3 <- psi*(theta.best[2]^2 *(Dt)^2/(theta.best[2]^2 *(Dt)^2+theta.best[4])/theta.best[4]+sqrt(theta.best[2]^2 *(Dt)^2+1)*theta.best[3]*Ds/2/(theta.best[2]^2 *(Dt)^2+theta.best[4])^(3/2))
    
    J0theta = matrix(0,5,5)
    g2 = (Gaminv%*%psic)
    g3 = (Gaminv%*%psi1)
    g4 = (Gaminv%*%psi2)
    g5 = (Gaminv%*%psi3)
    J0theta[1,1] = P/2/sigma2^2
    J0theta[2,2] = sum(g2*t(g2))/2
    J0theta[3,3] = sum(g3*t(g3))/2
    J0theta[4,4] = sum(g4*t(g4))/2
    J0theta[5,5] = sum(g5*t(g5))/2
    J0theta[2,3] = J0theta[3,2] = sum(g2*t(g3))/2
    J0theta[2,4] = J0theta[4,2] = sum(g2*t(g4))/2
    J0theta[2,5] = J0theta[5,2] = sum(g2*t(g5))/2
    J0theta[3,4] = J0theta[4,3] = sum(g3*t(g4))/2
    J0theta[3,5] = J0theta[5,3] = sum(g3*t(g5))/2
    J0theta[4,5] = J0theta[5,4] = sum(g4*t(g5))/2
    J0theta[2,1] = J0theta[1,2] = sum(diag(1/sigma2,P)*t(g2))/2
    J0theta[3,1] = J0theta[1,3] = sum(diag(1/sigma2,P)*t(g3))/2
    J0theta[4,1] = J0theta[1,4] = sum(diag(1/sigma2,P)*t(g4))/2
    J0theta[5,1] = J0theta[1,5] = sum(diag(1/sigma2,P)*t(g5))/2
  }
  if(CovStructure %in% c(3)){
    psic <- - psi / (1-theta.best[1]) # to c
    diag(psic) = 0
    psia <- -3 * psi * (Dt)^2 * theta.best[2] / (theta.best[2]^2 *(Dt)^2 + 1) # to ct
    psia[which(Dt==0,arr.ind = TRUE)] = 0
    psib <- - psi * Ds # to cs
    psib[which(Ds==0,arr.ind = TRUE)] = 0
  }
  if(CovStructure %in% c(4:8)){
    psic <- diag(P)
    psia <- -3 * psi * (Dt)^2 * theta.best[2] / (theta.best[2]^2 *(Dt)^2 + 1) # to ct
    psia[which(Dt==0,arr.ind = TRUE)] = 0
    psib <- - psi * Ds # to cs
    psib[which(Ds==0,arr.ind = TRUE)] = 0
  }
  if(CovStructure == 3){
    J0theta = matrix(0,4,4)
    g1 = (Gaminv%*%psic)
    g2 = (Gaminv%*%psia)
    g3 = (Gaminv%*%psib)
    
    J0theta[1,1] = sum(g1*t(g1))/2
    J0theta[2,2] = sum(g2*t(g2))/2
    J0theta[3,3] = sum(g3*t(g3))/2
    J0theta[4,4] = P/2/sigma2^2
    
    J0theta[1,2] = J0theta[2,1] = sum(g1*t(g2))/2
    J0theta[1,3] = J0theta[3,1] = sum(g1*t(g3))/2
    J0theta[2,3] = J0theta[3,2] = sum(g2*t(g3))/2
    J0theta[1,4] = J0theta[4,1] = sum(diag(1/sigma2,P)*t(g1))/2
    J0theta[2,4] = J0theta[4,2] = sum(diag(1/sigma2,P)*t(g2))/2
    J0theta[3,4] = J0theta[4,3] = sum(diag(1/sigma2,P)*t(g3))/2
  }
  if(CovStructure == 4){
    tmp1 = matrix(0,P,P)
    for(i in 1:P){
      for(j in 1:P){
        tmp1[i,j] = 1/(theta.best[4] + 1/loctim[i, 3]) + 1/(theta.best[4] + 1/loctim[j, 3])
      }
    }
    psid = tmp1 * psi
    diag(psid) = 2*loctim[,3]*DD1
    
    J0theta = matrix(0,5,5)
    g1 = (Gaminv%*%psic)
    g2 = (Gaminv%*%psia)
    g3 = (Gaminv%*%psib)
    g4 = (Gaminv%*%psid)
    
    J0theta[1,1] = sum(g1*t(g1))/2
    J0theta[2,2] = sum(g2*t(g2))/2
    J0theta[3,3] = sum(g3*t(g3))/2
    J0theta[4,4] = sum(g4*t(g4))/2
    J0theta[5,5] = P/2/sigma2^2
    
    J0theta[1,2] = J0theta[2,1] = sum(g1*t(g2))/2
    J0theta[1,3] = J0theta[3,1] = sum(g1*t(g3))/2
    J0theta[1,4] = J0theta[4,1] = sum(g1*t(g4))/2
    J0theta[2,3] = J0theta[3,2] = sum(g2*t(g3))/2
    J0theta[2,4] = J0theta[4,2] = sum(g2*t(g4))/2
    J0theta[3,4] = J0theta[4,3] = sum(g3*t(g4))/2
    J0theta[1,5] = J0theta[5,1] = sum(diag(1/sigma2,P)*t(g1))/2
    J0theta[2,5] = J0theta[5,2] = sum(diag(1/sigma2,P)*t(g2))/2
    J0theta[3,5] = J0theta[5,3] = sum(diag(1/sigma2,P)*t(g3))/2
    J0theta[4,5] = J0theta[5,4] = sum(diag(1/sigma2,P)*t(g4))/2
  }
  if(CovStructure == 5){
    tmp1 = matrix(0,P,P)
    for(i in 1:P){
      for(j in 1:P){
        tmp1[i,j] = 1/(theta.best[4] + (theta.best[5]*loctim[i,1] + theta.best[6]*loctim[i,2] + 1)/loctim[i,3]) +
          1/(theta.best[4] + (theta.best[5]*loctim[j,1] + theta.best[6]*loctim[j,2] + 1)/loctim[j,3])
      }
    }
    psid = tmp1 * psi
    diag(psid) = 2*loctim[,3]*DD1
    
    tmp2 = matrix(0,P,P)
    for(i in 1:P){
      for(j in 1:P){
        tmp2[i,j] = 1/(theta.best[5] + (theta.best[4]*loctim[i,3] + theta.best[6]*loctim[i,2] + 1)/loctim[i,1]) +
          1/(theta.best[5] + (theta.best[4]*loctim[j,3] + theta.best[6]*loctim[j,2] + 1)/loctim[j,1])
      }
    }
    psie = tmp2 * psi
    diag(psie) = 2*loctim[,1]*DD1
    
    tmp3 = matrix(0,P,P)
    for(i in 1:P){
      for(j in 1:P){
        tmp3[i,j] = 1/(theta.best[6] + (theta.best[4]*loctim[i,3] + theta.best[5]*loctim[i,1] + 1)/loctim[i,2]) +
          1/(theta.best[6] + (theta.best[4]*loctim[j,3] + theta.best[5]*loctim[j,1] + 1)/loctim[j,2])
      }
    }
    psif = tmp3 * psi
    diag(psif) = 2*loctim[,2]*DD1
    
    J0theta = matrix(0,7,7)
    g1 = (Gaminv%*%psic)
    g2 = (Gaminv%*%psia)
    g3 = (Gaminv%*%psib)
    g4 = (Gaminv%*%psid)
    g5 = (Gaminv%*%psie)
    g6 = (Gaminv%*%psif)
    
    J0theta[1,1] = sum(g1*t(g1))/2
    J0theta[2,2] = sum(g2*t(g2))/2
    J0theta[3,3] = sum(g3*t(g3))/2
    J0theta[4,4] = sum(g4*t(g4))/2
    J0theta[5,5] = sum(g5*t(g5))/2
    J0theta[6,6] = sum(g6*t(g6))/2
    J0theta[7,7] = P/2/sigma2^2
    
    J0theta[1,2] = J0theta[2,1] = sum(g1*t(g2))/2
    J0theta[1,3] = J0theta[3,1] = sum(g1*t(g3))/2
    J0theta[1,4] = J0theta[4,1] = sum(g1*t(g4))/2
    J0theta[1,5] = J0theta[5,1] = sum(g1*t(g5))/2
    J0theta[1,6] = J0theta[6,1] = sum(g1*t(g6))/2
    
    J0theta[2,3] = J0theta[3,2] = sum(g2*t(g3))/2
    J0theta[2,4] = J0theta[4,2] = sum(g2*t(g4))/2
    J0theta[2,5] = J0theta[5,2] = sum(g2*t(g5))/2
    J0theta[2,6] = J0theta[6,2] = sum(g2*t(g6))/2
    
    J0theta[3,4] = J0theta[4,3] = sum(g3*t(g4))/2
    J0theta[3,5] = J0theta[5,3] = sum(g3*t(g5))/2
    J0theta[3,6] = J0theta[6,3] = sum(g3*t(g6))/2
    
    J0theta[4,5] = J0theta[5,4] = sum(g4*t(g5))/2
    J0theta[4,6] = J0theta[6,4] = sum(g4*t(g6))/2
    
    J0theta[5,6] = J0theta[6,5] = sum(g5*t(g6))/2
    
    J0theta[1,7] = J0theta[7,1] = sum(diag(1/sigma2,P)*t(g1))/2
    J0theta[2,7] = J0theta[7,2] = sum(diag(1/sigma2,P)*t(g2))/2
    J0theta[3,7] = J0theta[7,3] = sum(diag(1/sigma2,P)*t(g3))/2
    J0theta[4,7] = J0theta[7,4] = sum(diag(1/sigma2,P)*t(g4))/2
    J0theta[5,7] = J0theta[7,5] = sum(diag(1/sigma2,P)*t(g5))/2
    J0theta[6,7] = J0theta[7,6] = sum(diag(1/sigma2,P)*t(g6))/2
  }
  if(CovStructure == 6)  {
    tmp1 = matrix(0,P,P)
    for(i in 1:P){
      for(j in 1:P){
        tmp1[i,j] = 1/(theta.best[4] + (theta.best[5]*(ifelse((loctim[i,3]-truncate.t)>0,(loctim[i,3]-truncate.t),0)) + 1)/loctim[i,3]) +
          1/(theta.best[4] + (theta.best[5]*(ifelse((loctim[j,3]-truncate.t)>0,(loctim[j,3]-truncate.t),0)) + 1)/loctim[j,3])
      }
    }
    psid = tmp1 * psi
    diag(psid) = 2*loctim[,3]*DD1
    
    tmp2 = matrix(0,P,P)
    for(i in 1:P){
      for(j in 1:P){
        tmp2[i,j] = (ifelse((loctim[i,3]-truncate.t)>0,(loctim[i,3]-truncate.t),0))/(1+theta.best[4]*loctim[i,3]+theta.best[5]*(ifelse((loctim[i,3]-truncate.t)>0,(loctim[i,3]-truncate.t),0)))+
          (ifelse((loctim[j,3]-truncate.t)>0,(loctim[j,3]-truncate.t),0))/(1+theta.best[4]*loctim[j,3]+theta.best[5]*(ifelse((loctim[j,3]-truncate.t)>0,(loctim[j,3]-truncate.t),0)))
      }
    }
    psie = tmp2 * psi
    diag(psie) = 2*(ifelse((loctim[,3]-truncate.t)>0,(loctim[,3]-truncate.t),0))*DD1
    
    J0theta = matrix(0,6,6)
    g1 = (Gaminv%*%psic)
    g2 = (Gaminv%*%psia)
    g3 = (Gaminv%*%psib)
    g4 = (Gaminv%*%psid)
    g5 = (Gaminv%*%psie)
    
    J0theta[1,1] = sum(g1*t(g1))/2
    J0theta[2,2] = sum(g2*t(g2))/2
    J0theta[3,3] = sum(g3*t(g3))/2
    J0theta[4,4] = sum(g4*t(g4))/2
    J0theta[5,5] = sum(g5*t(g5))/2
    J0theta[6,6] = P/2/sigma2^2
    
    J0theta[1,2] = J0theta[2,1] = sum(g1*t(g2))/2
    J0theta[1,3] = J0theta[3,1] = sum(g1*t(g3))/2
    J0theta[1,4] = J0theta[4,1] = sum(g1*t(g4))/2
    J0theta[1,5] = J0theta[5,1] = sum(g1*t(g5))/2
    J0theta[1,6] = J0theta[6,1] = sum(diag(1/sigma2,P)*t(g1))/2
    
    J0theta[2,3] = J0theta[3,2] = sum(g2*t(g3))/2
    J0theta[2,4] = J0theta[4,2] = sum(g2*t(g4))/2
    J0theta[2,5] = J0theta[5,2] = sum(g2*t(g5))/2
    J0theta[2,6] = J0theta[6,2] = sum(diag(1/sigma2,P)*t(g2))/2
    
    J0theta[3,4] = J0theta[4,3] = sum(g3*t(g4))/2
    J0theta[3,5] = J0theta[5,3] = sum(g3*t(g5))/2
    J0theta[3,6] = J0theta[6,3] = sum(diag(1/sigma2,P)*t(g3))/2
    
    J0theta[4,5] = J0theta[5,4] = sum(g4*t(g5))/2
    J0theta[4,6] = J0theta[6,4] = sum(diag(1/sigma2,P)*t(g4))/2
    
    J0theta[5,6] = J0theta[6,5] = sum(diag(1/sigma2,P)*t(g5))/2
  }
  if(CovStructure == 7){
    tmp1 = matrix(0,P,P)
    for(i in 1:P){
      for(j in 1:P){
        tmp1[i,j] = loctim[i,3]/DD1[i] + loctim[j,3]/DD1[j]
      }
    }
    psid = tmp1 * psi
    diag(psid) = 2*loctim[,3]*DD1
    
    tmp2 = matrix(0,P,P)
    for(i in 1:P){
      for(j in 1:P){
        tmp2[i,j] = loctim[i,3]^2/DD1[i] + loctim[j,3]^2/DD1[j]
      }
    }
    
    psie = tmp2 * psi
    diag(psie) = 2*loctim[,3]^2*DD1
    
    tmp3 = matrix(0,P,P)
    for(i in 1:P){
      for(j in 1:P){
        tmp3[i,j] = (ifelse((loctim[i,3]-truncate.t)>0,(loctim[i,3]-truncate.t)^2,0))/DD1[i] + 
          (ifelse((loctim[j,3]-truncate.t)>0,(loctim[j,3]-truncate.t)^2,0))/DD1[j]
      }
    }
    
    psif = tmp3 * psi
    diag(psif) = 2*(ifelse((loctim[,3]-truncate.t)>0,(loctim[,3]-truncate.t)^2,0))*DD1
    
    J0theta = matrix(0,7,7)
    g1 = (Gaminv%*%psic)
    g2 = (Gaminv%*%psia)
    g3 = (Gaminv%*%psib)
    g4 = (Gaminv%*%psid)
    g5 = (Gaminv%*%psie)
    g6 = (Gaminv%*%psif)
    
    J0theta[1,1] = sum(g1*t(g1))/2
    J0theta[2,2] = sum(g2*t(g2))/2
    J0theta[3,3] = sum(g3*t(g3))/2
    J0theta[4,4] = sum(g4*t(g4))/2
    J0theta[5,5] = sum(g5*t(g5))/2
    J0theta[6,6] = sum(g6*t(g6))/2
    J0theta[7,7] = P/2/sigma2^2
    
    J0theta[1,2] = J0theta[2,1] = sum(g1*t(g2))/2
    J0theta[1,3] = J0theta[3,1] = sum(g1*t(g3))/2
    J0theta[1,4] = J0theta[4,1] = sum(g1*t(g4))/2
    J0theta[1,5] = J0theta[5,1] = sum(g1*t(g5))/2
    J0theta[1,6] = J0theta[6,1] = sum(g1*t(g6))/2
    
    J0theta[2,3] = J0theta[3,2] = sum(g2*t(g3))/2
    J0theta[2,4] = J0theta[4,2] = sum(g2*t(g4))/2
    J0theta[2,5] = J0theta[5,2] = sum(g2*t(g5))/2
    J0theta[2,6] = J0theta[6,2] = sum(g2*t(g6))/2
    
    J0theta[3,4] = J0theta[4,3] = sum(g3*t(g4))/2
    J0theta[3,5] = J0theta[5,3] = sum(g3*t(g5))/2
    J0theta[3,6] = J0theta[6,3] = sum(g3*t(g6))/2
    
    J0theta[4,5] = J0theta[5,4] = sum(g4*t(g5))/2
    J0theta[4,6] = J0theta[6,4] = sum(g4*t(g6))/2
    
    J0theta[5,6] = J0theta[6,5] = sum(g5*t(g6))/2
    
    J0theta[1,7] = J0theta[7,1] = sum(diag(1/sigma2,P)*t(g1))/2
    J0theta[2,7] = J0theta[7,2] = sum(diag(1/sigma2,P)*t(g2))/2
    J0theta[3,7] = J0theta[7,3] = sum(diag(1/sigma2,P)*t(g3))/2
    J0theta[4,7] = J0theta[7,4] = sum(diag(1/sigma2,P)*t(g4))/2
    J0theta[5,7] = J0theta[7,5] = sum(diag(1/sigma2,P)*t(g5))/2
    J0theta[6,7] = J0theta[7,6] = sum(diag(1/sigma2,P)*t(g6))/2
  }
  if(CovStructure == 8){
    tmp1 = matrix(0,P,P)
    for(i in 1:P){
      for(j in 1:P){
        tmp1[i,j] = loctim[i,3]/DD1[i] + loctim[j,3]/DD1[j]
      }
    }
    psid = tmp1 * psi
    diag(psid) = 2*loctim[,3]*DD1
    
    tmp2 = matrix(0,P,P)
    for(i in 1:P){
      for(j in 1:P){
        tmp2[i,j] = loctim[i,3]^2/DD1[i] + loctim[j,3]^2/DD1[j]
      }
    }
    psie = tmp2 * psi
    diag(psie) = 2*loctim[,3]^2*DD1
    
    tmp3 = matrix(0,P,P)
    for(i in 1:P){
      for(j in 1:P){
        tmp3[i,j] = loctim[i,3]^3/DD1[i] + loctim[j,3]^3/DD1[j]
      }
    }
    psif = tmp3 * psi
    diag(psif) = 2*loctim[,3]^3*DD1
    
    tmp4 = matrix(0,P,P)
    for(i in 1:P){
      for(j in 1:P){
        tmp4[i,j] = (ifelse((loctim[i,3]-truncate.t)>0,(loctim[i,3]-truncate.t)^3,0))/DD1[i] + 
          (ifelse((loctim[j,3]-truncate.t)>0,(loctim[j,3]-truncate.t)^3,0))/DD1[j]
      }
    }
    psig = tmp4 * psi
    diag(psig) = 2*(ifelse((loctim[,3]-truncate.t)>0,(loctim[,3]-truncate.t)^3,0))*DD1
    
    J0theta = matrix(0,8,8)
    g1 = (Gaminv%*%psic)
    g2 = (Gaminv%*%psia)
    g3 = (Gaminv%*%psib)
    g4 = (Gaminv%*%psid)
    g5 = (Gaminv%*%psie)
    g6 = (Gaminv%*%psif)
    g7 = (Gaminv%*%psig)    
    
    J0theta[1,1] = sum(g1*t(g1))/2
    J0theta[2,2] = sum(g2*t(g2))/2
    J0theta[3,3] = sum(g3*t(g3))/2
    J0theta[4,4] = sum(g4*t(g4))/2
    J0theta[5,5] = sum(g5*t(g5))/2
    J0theta[6,6] = sum(g6*t(g6))/2
    J0theta[7,7] = sum(g7*t(g7))/2
    J0theta[8,8] = P/2/sigma2^2
    
    J0theta[1,2] = J0theta[2,1] = sum(g1*t(g2))/2
    J0theta[1,3] = J0theta[3,1] = sum(g1*t(g3))/2
    J0theta[1,4] = J0theta[4,1] = sum(g1*t(g4))/2
    J0theta[1,5] = J0theta[5,1] = sum(g1*t(g5))/2
    J0theta[1,6] = J0theta[6,1] = sum(g1*t(g6))/2
    J0theta[1,7] = J0theta[7,1] = sum(g1*t(g7))/2
    
    J0theta[2,3] = J0theta[3,2] = sum(g2*t(g3))/2
    J0theta[2,4] = J0theta[4,2] = sum(g2*t(g4))/2
    J0theta[2,5] = J0theta[5,2] = sum(g2*t(g5))/2
    J0theta[2,6] = J0theta[6,2] = sum(g2*t(g6))/2
    J0theta[2,7] = J0theta[7,2] = sum(g2*t(g7))/2
    
    J0theta[3,4] = J0theta[4,3] = sum(g3*t(g4))/2
    J0theta[3,5] = J0theta[5,3] = sum(g3*t(g5))/2
    J0theta[3,6] = J0theta[6,3] = sum(g3*t(g6))/2
    J0theta[3,7] = J0theta[7,3] = sum(g3*t(g7))/2
    
    J0theta[4,5] = J0theta[5,4] = sum(g4*t(g5))/2
    J0theta[4,6] = J0theta[6,4] = sum(g4*t(g6))/2
    J0theta[4,7] = J0theta[7,4] = sum(g4*t(g7))/2
    
    J0theta[5,6] = J0theta[6,5] = sum(g5*t(g6))/2
    J0theta[5,7] = J0theta[7,5] = sum(g5*t(g7))/2
    
    J0theta[6,7] = J0theta[7,6] = sum(g6*t(g7))/2
    
    J0theta[1,8] = J0theta[8,1] = sum(diag(1/sigma2,P)*t(g1))/2
    J0theta[2,8] = J0theta[8,2] = sum(diag(1/sigma2,P)*t(g2))/2
    J0theta[3,8] = J0theta[8,3] = sum(diag(1/sigma2,P)*t(g3))/2
    J0theta[4,8] = J0theta[8,4] = sum(diag(1/sigma2,P)*t(g4))/2
    J0theta[5,8] = J0theta[8,5] = sum(diag(1/sigma2,P)*t(g5))/2
    J0theta[6,8] = J0theta[8,6] = sum(diag(1/sigma2,P)*t(g6))/2
    J0theta[7,8] = J0theta[8,7] = sum(diag(1/sigma2,P)*t(g7))/2
  }
  
  Sigma = sigma2*psi
  
  result = list()
  result$beta = beta.est
  result$theta = c(sigma2,theta.best)
  result$J0betainv0 = J0betainv0
  result$J0theta = J0theta
  result$Gaminv = Gaminv
  if(all == TRUE){
    result$Sigma = Sigma
  }

  return(result)
}

