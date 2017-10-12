#' Partial Linear Model
#' 
#' Find least squares regression coefficients using only a subset of covariates
#' 
#' @param A the "active set", i.e. indices of covariates that will have non-zero coefficients
#' @param X matrix of predictor variables
#' @param y vector of response variables
#' 

partialLM <- function(X, y, A=1:ncol(X)) {
  if (length(A)>nrow(X)) {
    stop("number of covariates > number of observations, can't compute least squares estimates")
  }
  XA <- X[,A]
  coeff <- solve(t(XA) %*% XA) %*% (t(XA) %*% y) #calculate betas for active covariates
  beta <- numeric(ncol(X))
  beta[A] <- coeff
  beta[-A] <- 0
  beta
}