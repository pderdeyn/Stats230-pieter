#' HMM simulation occasional dishonest casino HMM
#'
#' @param P transition probabilities
#' @param E emission probabilities
#' @param v initial distribution
#' @param n number of simulated points
#' @return simulated values
#' @export
hmm_casino_simulate <- function(P,E,v,n=100) {
  x<-rep(0,n+1)
  y<-rep(0,n+1)
  x[1]<-sample(x=c(1,2),1,prob=v)
  y[1]<-sample(x=seq(6),1,prob=E[x[1],])
  for (i in seq(n)) {
    x[i+1]<-sample(x=c(1,2),1,prob=P[x[i],])
    y[i+1]<-sample(x=seq(6),1,prob=E[x[i],])
  }
  return(list(x,y))
}
