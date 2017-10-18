#' Forward stagewise Selection
#' 
#' Implements the Forward Stagewise Selection algorithm
#' 
#' @param X matrix of predictor variables
#' @param y vector of response variables
#' @param eps step size
#' @param tol stop when mu is within tol of y
#' @param N maximum number of steps
#' @param standardise should the data be centred and normalised?
#' @param intercept include an intercept term (i.e. constant covariate)
#' 

stagewise <- function(X, y, eps, tol=eps, N=1000, standardise=TRUE, intercept=FALSE) {

  ## Edit the data according to options
  if (intercept){
    X <- cbind(rep(1,nrow(X)), X) #include an intercept term (i.e. a constant covariate)
  }
  if (standardise){
    X <- scale(X) #centre and normalise X
    y <- y - mean(y) #centre y
  }
  
  n <- nrow(X) #number of observations
  p <- ncol(X) #number of covariates
  
  ## check inputs are sensible
  if (length(y)!=n){
    stop("number of observations on response not equal to number of observations on predictors")
  }
  if (length(eps) > 1){
    eps <- eps[1]
    warning("argument eps has length >1: only using first element")
  }
  if (length(tol) > 1){
    tol <- tol[1]
    warning("argument tol has length >1: only using first element")
  }
  if (length(N) > 1){
    N <- N[1]
    warning("argument N has length >1: only using first element")
  }
  if (tol < eps){
    tol <- eps
    warning("can't achieve tolerance lower than step size: setting tol=eps")
  }
  
  ## Initialise variables
  mu <- rep(0,n) #initial estimate is all zeroes
  beta <- rep(0,p) #initial coeficients beta is all zeroes
  M <- mu #matrix containing mu from each step (in columns)
  B <- beta #matrix containing beta from each step (in columns)
  t <- 0 #vector containing L1 norm of beta at each step
  J <- NA
  j <- 0
  delta <-0
  count <- 1

  ##let's go:
  while (any(abs(y-mu)>tol)){
      if (count>N){
        warning("maximum number of steps reached: consider increasing tol or N")
        break
      }
      j.old <- j
      delta.old <- delta
    
      j <- which.max(abs(t(X) %*% (y-mu))) #index of covariate most correlated to residual
      delta <- eps * sign(sum(X[,j] * (y-mu))) #step size & direction
      beta[j] <- beta[j] + delta
      mu <- mu + delta * X[,j]
      
      ## stop beta oscillating once converged:
      if (j.old == j && delta.old!=delta){
        break
      }
      
      M <- c(M,mu)
      B <- c(B,beta)
      J <- c(J,j)
      t <- c(t,sum(abs(beta)))
      count <- count + 1
  }
  
  M <- matrix(M, nrow=n)
  B <- matrix(B, nrow=p)
  
  #return some stuff:
  out <- list(beta=B, resid = apply(M,2,function(x){y-x}), t=t, j=J, method="Stagewise")
  class(out) <- "lars"
  out
  
}