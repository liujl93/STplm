#' @export
gt = function(t,h,FixT,KN,X,method=c("qt","noqt","plugin")) {
  qtfun = qt(FixT)
  switch(method,
         qt = t(kt(t,h,FixT,KN))%*%X/length(FixT)/qtfun(t),
         noqt = t(kt(t,h,FixT,KN))%*%X/length(FixT),
         plugin = t(kt(t,h,FixT,KN))%*%X/sum(kt(t,h,FixT,KN)))
}