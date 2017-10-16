#' Lasso
#' 
#' Implements the Lasso algorithm
#' 
#' @param X matrix of predictor variables
#' @param y vector of response variables
#' @param t shrinkage parameter
#' 

lasso<-function(X,y,t){
  n <- nrow(X) #number of observations
  p <- ncol(X) #number of covariates
  if (length(y)!=n){stop("number of observations on response not equal to number of observations on predictors")}
  ind <- 1:p
  
  M <- matrix(NA, nrow=n, ncol=p+1) #matrix containing mu from each step (in columns)
  B <- matrix(NA, nrow=p, ncol=p+1) #matrix containing beta from each step (in columns)
  M[,1] <- rep(0,n) #initial estimate is all zeroes
  B[,1] <- rep(0,p) #initial coeficients beta is all zeroes
  mu <- rep(0,n)
  
  # We use the algorithm detailed in Tibshirani 1994 
  # That uses the Kuhn-Tucker conditions to sequentially find feasible solutions
}