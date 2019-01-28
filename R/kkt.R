#' @export
kkt = function(t,h,FixT,KN,psi=NULL,sigma2,method=c("qt","noqt","plugin")){
  qtfun = qt(FixT)
  if(is.null(psi)){
    # changed 04/19/2017
    ret = switch(method,
                 qt = t(kt(t,h,FixT,KN))%*%kt(t,h,FixT,KN)*sigma2/length(FixT)^2/qtfun(t)^2,
                 noqt = t(kt(t,h,FixT,KN))%*%kt(t,h,FixT,KN)*sigma2/length(FixT)^2,
                 plugin = t(kt(t,h,FixT,KN))%*%kt(t,h,FixT,KN)*sigma2/(sum(kt(t,h,FixT,KN)))^2)}
    #ret = t(kt(t,h,FixT,KN))%*%kt(t,h,FixT,KN)*sigma2/h^2/length(FixT)^2}
  else {
    ret = switch(method,
                 qt = t(kt(t,h,FixT,KN))%*%psi%*%kt(t,h,FixT,KN)*sigma2/length(FixT)^2/qtfun(t)^2,
                 noqt = t(kt(t,h,FixT,KN))%*%psi%*%kt(t,h,FixT,KN)*sigma2/length(FixT)^2,
                 plugin = t(kt(t,h,FixT,KN))%*%psi%*%kt(t,h,FixT,KN)*sigma2/(sum(kt(t,h,FixT,KN)))^2)
    }
  #ret = t(kt(t,h,FixT,KN))%*%psi%*%kt(t,h,FixT,KN)*sigma2/h^2/length(FixT)^2}
  return(ret)
}

kkt1 = function(t,h,FixT,KN,psi=NULL,sigma2,method=c("qt","noqt","plugin")){
  qtfun = qt(FixT)
  if(is.null(psi)){
    # changed 04/19/2017
    ret = switch(method,
                 qt = t(kt1(t,h,FixT,KN))%*%kt1(t,h,FixT,KN)*sigma2/length(FixT)^2/qtfun(t)^2,
                 noqt = t(kt1(t,h,FixT,KN))%*%kt1(t,h,FixT,KN)*sigma2/length(FixT)^2,
                 plugin = t(kt1(t,h,FixT,KN))%*%kt1(t,h,FixT,KN)*sigma2/(sum(kt(t,h,FixT,KN)))^2)}
  #ret = t(kt(t,h,FixT,KN))%*%kt(t,h,FixT,KN)*sigma2/h^2/length(FixT)^2}
  else {
    ret = switch(method,
                 qt = t(kt1(t,h,FixT,KN))%*%psi%*%kt1(t,h,FixT,KN)*sigma2/length(FixT)^2/qtfun(t)^2,
                 noqt = t(kt1(t,h,FixT,KN))%*%psi%*%kt1(t,h,FixT,KN)*sigma2/length(FixT)^2,
                 plugin = t(kt1(t,h,FixT,KN))%*%psi%*%kt1(t,h,FixT,KN)*sigma2/(sum(kt(t,h,FixT,KN)))^2)
  }
  #ret = t(kt(t,h,FixT,KN))%*%psi%*%kt(t,h,FixT,KN)*sigma2/h^2/length(FixT)^2}
  ret = switch(KN,
               GK = ret,
               K2 = ret/2.25)
  return(ret)
}