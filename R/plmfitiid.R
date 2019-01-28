#' @import compiler
#' @export
PLM.fit.iid <- function(training, Type = 2, h0 = 0.005, h1 = 0.1, nlam = 30,beta.true=NULL) {
  loctim = as.matrix(training$loctim)
  P = nrow(loctim)
  X = training$X
  Z = training$Z

  if (Type == 1) { KERN1 = "GK" ;KERN2 = "GK"
  }else if (Type == 2) { KERN1 = "K2" ;KERN2 = "K2"
  }else if (Type == 3) { KERN1 = "K2" ;KERN2 = "GK"
  }else if (Type == 4) { KERN1 = "GK" ;KERN2 = "K2"
  }else stop("Type Not Defined.")


  bw = exp(seq(log(h0), log(h1), length = nlam))
  cv = gcv = numeric(nlam)

  if(is.null(beta.true)){
    result.ori = spl.iid(training, h0, KERN1)
    beta.ori = result.ori$beta.hat
  }else{
    beta.ori = beta.true
  }

  for (i in 1:nlam) {
    bw.fit = bw.cv.iid(training, bw[i], KERN2, beta.ori)
    cv[i] = bw.fit$CV
    gcv[i] = bw.fit$GCV
  }

  bw1 = bw[which.min(cv)]
  fit1 = spl.iid(training, bw1, KERN2)
  fit1$bandwidth = bw1

  bw2 = bw[which.min(gcv)]
  if(bw2 == bw1){
    fit2=fit1
  }else{
    fit2 = spl.iid(training, bw2, KERN2)
    fit2$bandwidth = bw2
  }
  ret = list()
  ret$cv = fit1
  ret$gcv = fit2


  return(ret)

}
