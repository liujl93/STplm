#' @export
qt = function(t){
  den = bkde(t,range.x = c(0,1))
  qtfun = approxfun(den$x, den$y, rule = 2:2)
  return(qtfun)
}