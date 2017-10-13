#' Partial Linear Model
#' 
#' Find least squares regression coefficients using only a subset of covariates
#'
#' @param X matrix of predictor variables
#' @param y vector of response variables
#' @param A the "active set", i.e. indices of covariates that will have non-zero coefficients
#' 

partialLM <- function(X, y, A=1:ncol(X)) {
  if (length(A)>nrow(X)) {
    stop("number of covariates > number of observations, can't compute least squares estimates")
  }
  X.A <- X[,A]
  coeff <- solve(t(X.A) %*% X.A) %*% (t(X.A) %*% y) #calculate betas for active covariates
  
  beta <- numeric(ncol(X))
  beta[A] <- coeff #put active coefficients in correct places
  beta[-A] <- 0 #fill inactive betas with zeroes
  
  beta
}