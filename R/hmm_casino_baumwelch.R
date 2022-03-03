#' HMM parameter estimation using forward and backward subsitution
#'
#' @param E emission probabilities for fair die
#' @param x hidden states
#' @param y observations
#' @param iter max iterations
#' @return all unkown parameters
#' @export
hmm_casino_baumwelch <- function(e,x,y,iter=1000,eps=0.05) {
  E <- matrix(c(e,rep(1/6,6)),nrow=2,byrow=TRUE)
  v <- c(0.5,0.5)
  P <- matrix(c(0.5,0.5,0.5,0.5),nrow=2)
  n <- length(y)
  sx <- dim(E)[1]
  sy <- dim(E)[2]
  lls <- rep(0,iter)
  for (k in seq(iter)) {

    # get marginals and forward/backward probs for updates

    fb_result<-hmm_casino_forwardbackward(P,E,v,y)
    gamma<-fb_result[[3]] # marginals
    a<-fb_result[[1]]
    b<-fb_result[[2]]

    # get log likelihood for previous params
    for (t in seq(n)) {
      lls[k]<-lls[k]+log(a[,t]%*%b[,t])
    }

    if (k>1) {
      if (abs(lls[k]-lls[k-1])<eps) {
        break
      }
    }

    # E update

    for (m in seq(sy)) {
      E[2,m]<-(gamma[2,]%*%(y==m))/sum(gamma[2,])
    }

    # P update

    g<-array(rep(0,sx*sx*n),c(sx,sx,n))
    for (i in seq(sx)) {
      for (j in seq(sx)) {
        for (t in seq(2,n)) {
          g[i,j,t]<-b[j,t]*E[j,y[t]]*P[i,j]*a[i,t-1]
        }
      }
    }
    for (t in seq(2,n)) {
      g[,,t]<-g[,,t]/sum(g[,,t])
    }

    for (i in seq(sx)) {
      for (j in seq(sx)) {
        P[i,j]<-sum(g[i,j,])
      }
    }
    for (j in seq(sx)) {
      P[,j]<-P[,j]/sum(P[,j])
    }

    # V update
    v<- gamma[,1]
  }
  return(list(lls[1:k],P,E,v))
}
