#' @import KernSmooth
#' @export
predict.plm.n <- function(object,training,test,Type,cross_section = FALSE,f.true=NULL,f.true0=NULL,DDnew=1,DDnew1=1,truncate.t=1){
  beta.hat = object$beta.hat
  f.hat = object$f.hat
  loctim = as.matrix(training$loctim)
  bandwidth = object$bandwidth
  CovStructure = object$CovStructure
  theta.hat = object$theta.hat
  Gaminv = object$Gaminv
    
  # Training Data
  P = nrow(loctim)
  X = as.matrix(training$X)
  Z = as.numeric(training$Z)
  Ds = as.matrix(training$Dloc)
  Dt = as.matrix(training$Dtim)
  # Test data
  sp0 = as.matrix(test$loctim[,-ncol(test$loctim)])
  time0 = as.numeric(test$loctim[, ncol(test$loctim)])
  P0 = nrow(sp0)
  
  Dloc0 = proxy::dist(sp0, loctim[, -ncol(loctim)]) * l
  Dtim0 = proxy::dist(time0, loctim[, ncol(loctim)]) * constt
  
  # Smooth Matrix
  if (Type == 1) { KERN1 = "GK" ;KERN2 = "GK"
  }else if (Type == 2) { KERN1 = "K2" ;KERN2 = "K2"
  }else if (Type == 3) { KERN1 = "K2" ;KERN2 = "GK"
  }else if (Type == 4) { KERN1 = "GK" ;KERN2 = "K2"
  }else stop("Type Not Defined.")
  theta.best = theta.hat[-1]
  if (CovStructure == 1) {
    DD1 = 1
    psi = (1-theta.best[1])*exp(-Ds/theta.best[2]-Dt/theta.best[3])
    diag(psi) = 1 
  }
  if (CovStructure == 2) {
    DD1 = 1
    psi = theta.best[4]/(theta.best[2]^2*Dt^2+1)^(1/2)/(theta.best[2]^2*Dt^2+theta.best[4])*exp(-theta.best[3]*Ds*((theta.best[2]^2*Dt^2+1)/(theta.best[2]^2*Dt^2+theta.best[4]))^(1/2))
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
  
  tmp = sqrt(DD1^2+ theta.hat[1])
  if(CovStructure %in% c(1:3)){ psi = psi/tmp^2
  }else{psi = psi/(tmp%*%t(tmp))}
  # ----------------------------------------------------- #    
  #    Inverse correlation matrix (Gamma/sigma^2)^(-1)    #
  # ----------------------------------------------------- #  
  U <- base::chol(psi)
  U.inv = backsolve(U,diag(1,nrow=P))
  Gaminv = U.inv%*%t(U.inv)
  
  if(is.null(f.true)){
    if(cross_section==FALSE){
      Smat0 = matrix(0,P0,P)
      for(i in 1:P0) Smat0[i,] = St(test$loctim[i,3],bandwidth,loctim[,3],KERN2)[1,]
      fhat0 = Smat0%*%(Z - X %*% beta.hat)
    } else{
      smat0 = St(test$loctim[1,3],bandwidth,loctim[,3],KERN2)[1,]
      fhat0 = as.numeric(smat0%*%(Z - X %*% beta.hat))
    }
  }
  
  if(is.null(object$CovStructure)){
    pred = test$X %*% beta.hat + fhat0
  }  else{
    if (CovStructure == 1) {
      DD1 = DD2 = 1
      Gam0 = (1 - theta.hat[2]) * exp(-Dloc0/theta.hat[3] - Dtim0/theta.hat[4])
    }
    if (CovStructure == 2) {
      DD1 = DD2 = 1
      Gam0 = (1-theta.hat[2])/(theta.hat[3]*(Dtim0)^2+1)*exp(-theta.hat[4]*Dloc0^2/(theta.hat[3]*(Dtim0)^2 + 1))
    }
    if (CovStructure == 3) {
      DD1 = DD2 = 1
      Gam0 <- 1/(theta.hat[3]^2*(Dtim0)^2+1)^(3/2)*exp(-theta.hat[4]*Dloc0)
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
      DD1 =1 + theta.hat[5]*time0 + theta.hat[6]*(ifelse((time0-truncate.t)>0,(time0-truncate.t),0))
      DD2 =1 + theta.hat[5]*(loctim[,3]) + theta.hat[6]*(ifelse((loctim[,3]-truncate.t)>0,(loctim[,3]-truncate.t),0))
      Gam0 <- DD1 %*% t(DD2)/(theta.hat[3]^2*(Dtim0)^2 + 1)^(3/2) * exp(-theta.hat[4] * Dloc0)
    }
    if (CovStructure == 7) {
      DD1 =1 + theta.hat[5]*time0 + theta.hat[6]*time0^2 + theta.hat[7]*(ifelse((time0-truncate.t)>0,(time0-truncate.t)^2,0))
      DD2 =1 + theta.hat[5]*(loctim[,3]) + theta.hat[6]*(loctim[,3])^2 + theta.hat[7]*(ifelse((loctim[,3]-truncate.t)>0,(loctim[,3]-truncate.t)^2,0))
      Gam0 <- DD1 %*% t(DD2)/(theta.hat[3]^2*(Dtim0)^2 + 1)^(3/2) * exp(-theta.hat[4] * Dloc0)
    }
    if (CovStructure == 8) {
      DD1 =1 + theta.hat[5]*time0 + theta.hat[6]*time0^2 + theta.hat[7]*time0^3 + theta.hat[8]*(ifelse((time0-truncate.t)>0,(time0-truncate.t)^3,0))
      DD2 =1 + theta.hat[5]*(loctim[,3]) + theta.hat[6]*(loctim[,3])^2 + theta.hat[7]*(loctim[,3])^3 + theta.hat[8]*(ifelse((loctim[,3]-truncate.t)>0,(loctim[,3]-truncate.t)^3,0))
      Gam0 <- DD1 %*% t(DD2)/(theta.hat[3]^2*(Dtim0)^2 + 1)^(3/2) * exp(-theta.hat[4] * Dloc0)
    }
    if (CovStructure == 9) {
      DD1 =1 + theta.hat[5]*(ifelse((time0-truncate.t)>0,1,0))
      DD2 =1 + theta.hat[5]*(ifelse((loctim[,3]-truncate.t)>0,1,0))
      Gam0 <- DD1 %*% t(DD2)/(theta.hat[3]^2*(Dtim0)^2 + 1)^(3/2) * exp(-theta.hat[4] * Dloc0)
    }
    if (CovStructure == 10) {
      DD1 =1 + theta.hat[5]*(ifelse((time0-truncate.t[1])>0,(time0-truncate.t[1]),0)) + theta.hat[6]*(ifelse((time0-truncate.t[2])>0,(time0-truncate.t[2]),0))
      DD2 =1 + theta.hat[5]*(ifelse((loctim[,3]-truncate.t[1])>0,(loctim[,3]-truncate.t[1]),0)) + theta.hat[6]*(ifelse((loctim[,3]-truncate.t[2])>0,(loctim[,3]-truncate.t[2]),0))
      Gam0 <- DD1 %*% t(DD2)/(theta.hat[3]^2*(Dtim0)^2 + 1)^(3/2) * exp(-theta.hat[4] * Dloc0)
    }
    if (CovStructure == 11) {
      DD1 = DDnew1    
      DD2 = DDnew
      Gam0 <- DD1 %*% t(DD2)/(theta.hat[3]^2*(Dtim0)^2 + 1)^(3/2) * exp(-theta.hat[4] * Dloc0)
    }
      tmp1 = sqrt(DD1^2+ theta.hat[1])
      tmp2 = sqrt(DD2^2+ theta.hat[1])
      if(CovStructure %in% c(1:3)){ Gam0 = Gam0/tmp1^2
      }else{Gam0 = Gam0/(tmp1%*%t(tmp2))}
    
    if(is.null(f.true)){
      pred = test$X %*% beta.hat + fhat0 + Gam0 %*% Gaminv %*% (Z - X %*% beta.hat - f.hat)
    }else{
      pred = test$X %*% beta.hat + f.true0 + Gam0 %*% Gaminv %*% (Z - X %*% beta.hat - f.true)
    }
  }
  return(pred)
}

