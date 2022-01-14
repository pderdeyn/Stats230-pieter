#' Multiple two matrices and a vector
#'
#' @param A a square matrix
#' @param B a square matrix with the same dim as A
#' @param x a vector with length equal to length of A
#' @param left a boolean option that determines which product is computed first
#' @examples
#' mat_mult(cbind(c(1,2),c(3,4)),rbind(c(1,2),(3,4)),c(1,2),True)
#' @return matrix product \eqn{A \times B \times x}
#' @export
mat_mult <- function(A,B,x,left=TRUE){
  prod <- 1
  if (left) {
    prod <- A %*% B
    prod <- prod %*% x
  }
  else {
    prod <- B %*% x
    prod <- A %*% prod
  }
  return(prod)
}
