#' Simulation of MVN random vectors
#'
#' @param mu n-dim vector, mean
#' @param sigma nxn positive definite matrix, covariance matrix
#' @param N number of relizations from MVN(mu,sigma)
#' @examples
#' mat_mult(cbind(c(1,2),c(3,4)),rbind(c(1,2),(3,4)),c(1,2),True)
#' @return N realizations from MVN distribution}
#' @export
mvn_chol <- function(mu,sigma,N){
  L<-chol(sigma)
  k<-length(mu)
  sim <- matrix(data=rnorm(N*k), nrow=k, ncol=N)
  sim <- mu + L%*%sim
  return(sim)
}
