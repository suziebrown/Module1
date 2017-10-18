#' Forward Selection 
#' 
#' Implements the Classic Forward Selection algorithm
#' 
#' @param X matrix of predictor variables
#' @param y vector of response variables
#' @param standardise should the data be centred and normalised?
#' @param intercept include an intercept term (i.e. constant covariate)
#' 
#' @return an object of class "lars"
#'
#' @details Classic forward selection is a naive model selection tool for linear regression. It tends to be overly greedy and produce model trees that are very susceptible to perturbations of the data. It is not recommended to use this method for anything serious, but it might be fun as a toy.
#' 
#' @export forward
#' 

forward <- function(X, y, standardise=T, intercept=F) {
  
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
  
  ## Check inputs are sensible
  if (length(y)!=n){
    stop("number of observations on response not equal to number of observations on predictors")
  }
  if (p>n) { #LSE can't be calculated with >n covariates, return NAs
    warning("model parameters can only be estimated with up to n covariates, only returning first n models")
  }

  ## Initialise variables
  A <- numeric(0) #active set
  M <- matrix(NA, nrow=n, ncol=p+1) #matrix containing mu from each step (in columns)
  B <- matrix(NA, nrow=p, ncol=p+1) #matrix containing beta from each step (in columns)
  M[,1] <- rep(0,n) #initial estimate is all zeroes
  B[,1] <- rep(0,p) #initial coeficients beta is all zeroes
  t <- 0 #vector containing L1 norm of beta at each step
  mu <- rep(0,n)
  J <- NA
  ind <- 1:p

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
    t[i] <- sum(abs(beta))
    J <- c(J,j)
  }
  
  out <- list(beta=B, resid = apply(M,2,function(x){y-x}), t=t, j=J, method="Forward")
  class(out) <- "lars"
  out
}