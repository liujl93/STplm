#' @import compiler
#' @export
bw.cv.iid <- function(training, h, KERN, beta.fixed) 
{
  loctim = as.matrix(training$loctim)
  P = nrow(loctim)
  X = training$X
  Z = training$Z
  Smat = matrix(0, P, P)
  for (i in 1:P) {
      Smat[i, ] = St(loctim[i, 3], h, loctim[, 3], KERN)[1,]
  }
  Zprime = (diag(P) - Smat) %*% Z
  Xprime = (diag(P) - Smat) %*% X
  eprime = Zprime - Xprime %*% beta.fixed

  CV = mean((eprime/(1 - diag(Smat)))^2)
  GCV = mean((eprime/(1 - 1/P * sum(diag(Smat))))^2)

  result = list()
  result$CV = CV
  result$GCV = GCV
  return(result)
}
