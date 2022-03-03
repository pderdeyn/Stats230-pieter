#' fit logistic regression model to data
#'
#' @param X data
#' @param y labels
#' @param method optimization method, by default gradient descent
#' @param iter max number of iterations
#' @param stepsize gradient descent stepsize
#' @param tol stopping criteria for log likelihoods
#' @return m model coefficients,
#' @export
logreg_fit <- function(X,y,method='gradient descent',iter=1000,stepsize=0.001,tol=0.01){
  if (method=='gradient descent') {
    ndim <- dim(X)[2]
    theta<-numeric(ndim)
    theta_history<-matrix(rep(0,ndim*iter),nrow=ndim)
    dLL <- numeric(ndim)
    LL <- numeric(iter+1)

    for (row in 1:nrow(X)) {
      LL[1] <- LL[1] + log_likelihood(theta,X[row,],y[row])
    }
    for (i in seq(iter)) {
      dLL<-0
      for (row in 1:nrow(X)) {
        dLL<-dLL+gradient(theta,X[row,],y[row])
      }
      theta<-theta+stepsize*dLL
      theta_history[,i]<-theta

      for (row in 1:nrow(X)) {
        LL[i+1]<-LL[i+1]+log_likelihood(theta,X[row,],y[row])
      }
      if (abs(LL[i+1]-LL[i]) < tol) {
        break
      }
    }
  }

  else if (method=='newton-raphson') {
    ndim <- dim(X)[2]
    theta<-numeric(ndim)
    theta_history<-matrix(rep(0,ndim*iter),nrow=ndim)
    dLL <- numeric(ndim)
    LL <- numeric(iter+1)

    for (row in 1:nrow(X)) {
      LL[1] <- LL[1] + log_likelihood(theta,X[row,],y[row])
    }
    H<-matrix(rep(0,length(theta)^2),nrow=length(theta))
    for (i in seq(iter)) {
      dLL <- numeric(ndim)
      for (row in 1:nrow(X)) {
        H <- H + hessian(theta,X[row,],y[row])
        dLL<-dLL+gradient(theta,X[row,],y[row])
      }
      Hinv <- solve(H)
      theta <- theta + Hinv %*% dLL
      theta_history[,i]<-theta

      for (row in 1:nrow(X)) {
        LL[i+1]<-LL[i+1]+log_likelihood(theta,X[row,],y[row])
      }
      if (abs(LL[i+1]-LL[i]) < tol) {
        break
      }
    }

  }
  theta_confidence <- rep(0,ndim)
  theta_history<-theta_history[,(i-10):i]
  for (f in seq(ndim)) {
    theta_confidence[f] <- qnorm(0.975)*sd(theta_history[f,])/sqrt(11)
  }
  return(list(theta,LL[1:(i+1)],theta_confidence))
}

sigmoid <- function(theta,x) {
  return(1/(1+exp(-sum(theta*x))))
}

log_likelihood <- function(theta,x,y) {
  sig<-sigmoid(theta,x)
  return(y*log(sig)+(1-y)*log(1-sig))
}

gradient <- function(theta,x,y) {
  sig<-sigmoid(theta,x)
  return(x*rep(y-sig,length(y)))
}

hessian <- function(theta,x,y) {
  sig<-sigmoid(theta,x)
  H<-matrix(rep(0,length(theta)^2),nrow=length(theta))
  H<- sig*(1-sig)*(x%*%t(x))
}
