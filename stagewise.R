#' Forward stagewise Selection
#' 
#' Implements the Forward Stagewise Selection algorithm
#' 
#' @param X matrix of predictor variables
#' @param y vector of response variables
#' @param eps step size
#' @param N number of steps
#' 

stagewise <- function(X, y, eps, N=5) {
  n <- nrow(X) #number of observations
  p <- ncol(X) #number of covariates
  if (length(y)!=n){stop("number of observations on response not equal to number of observations on predictors")}

  M <- matrix(NA, nrow=n, ncol=N) #matrix containing mu from each step (in columns)
  B <- matrix(NA, nrow=p, ncol=N) #matrix containing beta from each step (in columns)
  mu <- rep(0,n) #initial estimate is all zeroes
  beta <- rep(0,p) #initial coeficients beta is all zeroes
  M[,1] <- mu
  B[,1] <-beta
  js <- numeric(N)
  js[1] <-NA

  for (i in 2:N) {
      j <- which.max(t(X) %*% (y-mu)) #index of covariate most correlated to residual
      delta <- eps * sign(t(X[,j]) %*% (y-mu)) #step size & direction
      beta[j] <- beta[j] + delta
      mu <- mu + delta * X[,j]
      M[,i] <- mu
      B[,i] <- beta
      js[i] <- j
  }
  list(coeffs=B, predicts=M, moved=js)
}