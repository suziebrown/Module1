#' Lasso wrapper
#' 
#' Runs the Lasso algorithm for a range of t
#' And outputs the results in the correct form
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
  p <- ncol(X) #number of covariates
  # Initialize matrices for the betas and the L1 norm t
  beta_mat<-rep(0,p)
  t_start<-0
  t_vec<-t_start
  for (t1 in t1_range){
    beta_new<-lasso_step(X, y, t1, eps, N, standardise, intercept)
    beta_mat<-cbind(beta_mat,beta_new)
    t_vec<-c(t_vec,sum(abs(beta_new)))
  }
  beta_mat<-matrix(beta_mat,p,length(t1_range)+1)
  
  output<-list(beta=beta_mat, resid=NA, t=t_vec, j=NA, method="Lasso")
  class(output)<-"lars"
  output
}
