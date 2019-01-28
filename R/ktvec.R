#' @export
kt = function(t,h,FixT,KN) {
  switch(KN,
         GK = c(dnorm((FixT-t)/h)/h),
         K2 = c(K2((FixT-t)/h))/h)
}

kt1 = function(t,h,FixT,KN) {
  switch(KN,
         GK = c((FixT-t)/h*dnorm((FixT-t)/h)/h),
         K2 = c((FixT-t)/h*K2((FixT-t)/h)/h))
}

kt2 <- function (t, h, FixT, KN)
{
  switch(KN, GK = c((FixT - t)/h * dnorm((FixT - t)/h)/h),
         K2 = c((FixT - t)/h * K2((FixT - t)/h))/h)
}
