#' Forward Selection 
#' 
#' Implements the Forward Selection algorithm
#' 
#' @param X matrix of predictor variables
#' @param y vector of response variables
#' 

forward <- function(X, y) {
  n <- nrow(X) #number of observations
  p <- ncol(X) #number of covariates
  ind <- 1:p

  A <- numeric(0) #active set
  M <- matrix(NA, nrow=p, ncol=p) #matrix containing mu at each step
  mu <- rep(0,p) #initial estimate is all zeroes
  
  if (p>n) { #LSE can't be calculated with >n covariates, return NAs
    M[-(1:n),] <- NA 
    warning("model parameters can only be estimated with up to n covariates, returning NA for higher-dimensional models")
  }

  for (i in seq_along(min(n,p))) {
    ind <- ind[-A] #inactive indices
    M[i,] <- mu 
    s <- sign(t(X[,-A]) %*% (y-mu))
    XA <- t(s*t(X[,-A]))
    j <- which.max(t(XA) %*% (y-mu)) #most correlated covariate not in active set
    j <- ind[j] #convert index with A-subset to true index
    A <- c(A,j) #add j to active set
    beta <- partialLM(X,y,A) #find regression coefficients on active set
    mu <- X %*% beta
  }
  M
}