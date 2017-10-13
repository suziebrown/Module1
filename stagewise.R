#' Forward stagewise Selection
#' 
#' Implements the Forward Stagewise Selection algorithm
#' 
#' @param X matrix of predictor variables
#' @param y vector of response variables
#' @param eps step size
#' @param tol stop when mu is within tol of y
#' @param N maximum number of steps
#' 

stagewise <- function(X, y, eps, tol=10*eps, N=1000) {
  
  #X <- rbind(X, rep(1, ncol(X))) #to include an intercept term(?)
  n <- nrow(X) #number of observations
  p <- ncol(X) #number of covariates
  
  #check inputs are sensible:
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
  if (tol < 10*eps){
    tol <- 10*eps
    warning("can't achieve tolerance of same magnitude as step size: setting tol=10*eps")
  }
  
  #initialise variables:
  mu <- rep(0,n) #initial estimate is all zeroes
  beta <- rep(0,p) #initial coeficients beta is all zeroes
  M <- mu #matrix containing mu from each step (in columns)
  B <- beta #matrix containing beta from each step (in columns)
  count <- 1

  #let's go:
  while (any(abs((y-mu))>tol)){
      if (count>N){
        warning("maximum number of steps reached: consider increasing tol or increasing N")
        break
      }
      j <- which.max(abs(t(X) %*% (y-mu))) #index of covariate most correlated to residual
      delta <- eps * sign(sum(X[,j] * (y-mu))) #step size & direction
      beta[j] <- beta[j] + delta
      mu <- mu + delta * X[,j]
      M <- c(M,mu)
      B <- c(B,beta)
      count <- count + 1
  }
  
  M <- matrix(M, nrow=n)
  B <- matrix(B, nrow=p)
  #return some stuff:
  list(coeffs=B, predicts=M)
  
}