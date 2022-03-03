#' HMM parameter estimation using forward and backward subsitution
#'
#' @param P transition probabilities
#' @param E emission probabilities
#' @param v initial distribution
#' @param y observations
#' @return forward and backward probs, and hidden state marginal probs
#' @export
hmm_casino_forwardbackward <- function(P,E,v,y) {
  n<-length(y)
  s<-dim(E)[1]
  b<-matrix(rep(0,n*s),ncol=n)
  b[,n]<-rep(1,s)
  for (t in seq(n-1,2)) {
    for (i in seq(s)) {
      for (j in seq(s)) {
        b[i,t]=b[i,t]+P[j,i]*E[j,y[t+1]]*b[j,t+1]
      }
    }
  }
  for (i in seq(s)) {
    for (j in seq(s)) {
      b[i,1]=b[i,1]+v[j]*E[j,y[2]]*b[j,2]
    }
  }

  a<-matrix(rep(0,n*s),ncol=n)
  a[,1] <- v
  for (i in seq(s)) {
    a[i,1] <- a[i,1] * E[i,y[1]]
  }

  for (t in seq(n-1)) {
    for (i in seq(s)) {
      a[i,t+1] <- E[i,y[t+1]]*(a[,t] %*% P[,i])
    }
  }
  m<-a*b
  for (i in seq(n)) {
    m[,i]<-m[,i]/sum(m[,i])
  }

  return(list(a,b,m))
}
