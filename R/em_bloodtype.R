#' EM algorithm for ABO blood type
#'
#' @param nA fraction of A blood type
#' @param nAB fraction of AB blood type
#' @param nB fraction of B blood type
#' @param nO fraction of O blood type
#' @param PA initial guess for A allele frequency
#' @param pB initial guess for B allele frequency
#' @param pO initial guess for O allele frequency
#' @param iter max number of iterations
#' @param epsilon stopping criteria for allele frequencies
#' @param delta stopping criteria for observed likelihood
#' @return population allele frequencies and observed data likelihood
#' @export
em_bloodtype <- function(nA,nAB,nB,nO,pA0,pB0,pO0,iter=1000,epsilon=0.01,delta=0.01) {
  mAA<-rep(0,iter)
  mAB<-nAB
  mAO<-rep(0,iter)
  mBB<-rep(0,iter)
  mBO<-rep(0,iter)
  mOO<-nO

  pA<-rep(0,iter)
  pB<-rep(0,iter)
  pO<-rep(0,iter)

  pA[1]<-pA0
  pB[1]<-pB0
  pO[1]<-pO0

  l_old <-(pA[1]^2+2*pA[11]*pO[1])^nA *
    (pB[1]^2+2*pB[1]*pO[1])^nB *
    (2*pA[1]*pB[1])^nAB*(pO[1]^2)^nO

  ll_old <- log(l_old)
  lls <- rep(0,iter+1)
  lls[1] <- ll_old

  for (i in seq(iter+1)) {
    # E-step
    mAA[i] <- nA*pA[i]^2/(pA[i]^2+2*pA[i]*pO[i])
    mAO[i] <- nA*2*pA[i]*pO[i]/(pA[i]^2+2*pA[i]*pO[i])
    mBB[i] <- nB*pB[i]^2/(pB[i]^2+2*pB[i]*pO[i])
    mBO[i] <- nB*2*pB[i]*pO[i]/(pB[i]^2+2*pB[i]*pO[i])

    # M-step
    pA[i+1]<-(2*mAA[i]+mAO[i]+mAB)/2
    pB[i+1]<-(2*mBB[i]+mBO[i]+mAB)/2
    pO[i+1]<-(2*mOO+mAO[i]+mBO[i])/2

    epsilon_current<-max(abs(pA[i+1]-pA[i]),
        abs(pB[i+1]-pB[i]),
        abs(pO[i+1]-pO[i]))

    l_new <- (pA[i+1]^2+2*pA[i+1]*pO[i+1])^nA *
      (pB[i+1]^2+2*pB[i+1]*pO[i+1])^nB *
      (2*pA[i+1]*pB[i+1])^nAB*(pO[i+1]^2)^nO
    ll_new <- log(l_new)
    delta_current<-abs(ll_new-ll_old)
    ll_old <- ll_new
    lls[i+1]<-ll_new

    if (i > 1000) {
      if (i %% 100 == 0) {
        print('eps')
        print(epsilon_current)
        print('delta')
        print(delta_current)
      }
    }

    if (epsilon_current < epsilon) {
      print('passed param tolerance')
      print(epsilon_current)
      if (delta_current < delta) {
        print('passed likelihood tolerance')
        print(delta_current)
        break
      }
    }
  }
  return(list(pA[1:(i+1)],pB[1:(i+1)],pO[1:(i+1)],l_new,lls[1:(i+1)]))
}
