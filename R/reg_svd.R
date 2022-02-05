#' Linear regression with SVD decomposition
#'
#' @param X covariates
#' @param y response
#' @return regresion coefficients
#' @export
reg_svd <- function(x,y){
  udv<-svd(x)
  u<-udv$u
  d<-udv$d
  d<-1/d
  D<-diag(d)
  v<-udv$v
  coef <- v %*% D %*% t(u) %*% y
  return(coef)
}
