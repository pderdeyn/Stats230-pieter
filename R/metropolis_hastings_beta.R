#' metropolis hastings for beta distribution wiht log normal porpo
#'
#' @param a parameter for beta
#' @param b parameter for beta
#' @param N length of markov chain
#' @param sigma for log normal proposal density
#' @return markov chain
#' @export
metropolis_hastings_beta <- function (a=4, b=2, N=1000, sigma=1)
{
  x<-rep(0,N)
  xcur<-rbeta(1,a,b)
  for (n in seq(N)) {
    #step <- rlnorm(1,meanlog=0,sigma)*sample(c(-1,1),1)
    #xprop<-xcur+step
    #acceptance<-dbeta(xprop,a,b)/dbeta(xcur,a,b)
    xprop <- rlnorm(1,meanlog=xcur,sigma)
    acceptance<-(dbeta(xprop,a,b)*dlnorm(xcur,xprop,sigma))/(dbeta(xcur,a,b)*dlnorm(xprop,xcur,sigma))
    u<-runif(1,0,1)
    if (u<=acceptance) {
      x[n]<-xprop
    }
    else {
      x[n]<-xcur
    }
    xcur<-x[n]

  }
  return(x)
}
