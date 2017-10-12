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
  count <- 0
  
  A <- numeric(0) #active set
  mu <- rep(0,p)

  for (count in seq_along(p)) {
    j <- which.max(X[-A]%*%(y-mu)) #most correlated covariate not in active set
    A <- c(A,j) #add j to active set
    beta <- partialLM(X,y,A) #find regression coefficients on active set
    # ... not finished yet
  }
}