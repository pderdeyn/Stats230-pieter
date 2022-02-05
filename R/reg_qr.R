#' Linear regression with QR decomposition
#'
#' @param X covariates
#' @param y response
#' @return regresion coefficients
#' @export
reg_qr <- function(x,y){
  QR <- qr(x)
  Q <- qr.Q(QR)
  R <- qr.R(QR)
  coef<-solve(R,t(Q)%*%y)
  return(coef)
}
