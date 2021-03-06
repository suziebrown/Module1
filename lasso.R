#' Lasso step
#' 
#' Implements the Lasso algorithm for a specific t
#' 
#' @param X matrix of predictor variables
#' @param y vector of response variables
#' @param t1 shrinkage parameter
#' @param eps cutoff for determining when a variable is identically 0
#' @param N maximum number of steps for numerical procedure
#' @param standardise option to centre and normalise the covariate matrix X
#' @param intercept option to add an intercept term to the covariate matrix
#' 


lasso<-function(X, y, t1_range, eps=1e-3, N=100, standardise=TRUE, intercept=FALSE){
  n <- nrow(X) #number of observations
  p <- ncol(X) #number of covariates
  # Initialize matrices for the mus, betas and the L1 norm t
  res_mat <- rep(0,n) # initialize matrix containing residuals from each step (in columns)
  beta_mat<-rep(0,p) # initialize matrix containing beta from each step (in columns)
  t_vec<-0 # vector containing t from each step
  for (t1 in t1_range){
    output<-lasso_step(X, y, t1, eps, N, standardise, intercept)
    
    beta_new<-output$beta
    res_new<-output$res
    beta_mat<-cbind(beta_mat,beta_new)
    res_mat<-cbind(res_mat,res_new)
    
    t_vec<-c(t_vec,sum(abs(beta_new)))
  }
  beta_mat<-matrix(beta_mat,p,length(t1_range)+1)
  
  
  output<-list(beta=beta_mat, resid=res_mat, t=t_vec, j=NA, method="Lasso")
  class(output)<-"lars"
  output
}
