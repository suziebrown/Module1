#' Forward Selection 
#' 
#' Implements the Forward Selection algorithm
#' 
#' @param X matrix of predictor variables
#' @param y vector of response variables
#' 

forward <- function(X, y) {
  X <- scale(X)
  n <- nrow(X) #number of observations
  p <- ncol(X) #number of covariates
  if (length(y)!=n){stop("number of observations on response not equal to number of observations on predictors")}
  ind <- 1:p

  A <- numeric(0) #active set
  M <- matrix(NA, nrow=n, ncol=p+1) #matrix containing mu from each step (in columns)
  B <- matrix(NA, nrow=p, ncol=p+1) #matrix containing beta from each step (in columns)
  M[,1] <- rep(0,n) #initial estimate is all zeroes
  B[,1] <- rep(0,p) #initial coeficients beta is all zeroes
  mu <- rep(0,n)
  
  if (p>n) { #LSE can't be calculated with >n covariates, return NAs
    warning("model parameters can only be estimated with up to n covariates, returning NA for higher-dimensional models")
  }

  for (i in 2:(min(n,p)+1)) {
    if(length(A)==0) { #create signed submatrix of X; case of empty A handled separately
      s <- sign(t(X) %*% (y-mu))
      X.notA <- rep(s, each=n) * X
    }
    else{
      s <- sign(t(X[,-A]) %*% (y-mu))
      X.notA <- rep(s, each=n) * X[,-A]
    }
    j <- which.max(t(X.notA) %*% (y-mu)) #most correlated covariate not in active set
    j <- ind[j] #convert index within active set to true index
    A <- c(A,j) #add j to active set
    ind <- setdiff(ind, A) #update inactive indices
    beta <- partialLM(X,y,A) #find regression coefficients on active set
    mu <- X %*% beta
    M[,i] <- mu 
    B[,i] <- beta
  }
  list(coeffs=B, predicts=M)
}