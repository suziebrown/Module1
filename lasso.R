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

lasso_step<-function(X, y, t1, eps=1e-3, N=100, standardise=TRUE, intercept=FALSE){
  n <- nrow(X) #number of observations
  p <- ncol(X) #number of covariates
  
  ## Edit the data according to options
  if (intercept){
    X <- cbind(rep(1,n), X) #include an intercept term (i.e. a constant covariate)
  }
  if (standardise){
    X <- scale(X) #centre and normalise X
    y <- y - mean(y) #centre y
  }
  
  ##check inputs are sensible:
  if (length(y)!=n){
    stop("number of observations on response not equal to number of observations on predictors")
  }
  if (length(eps) > 1){
    eps <- eps[1]
    warning("argument eps has length >1: only using first element")
  }
  if (length(N) > 1){
    N <- N[1]
    warning("argument N has length >1: only using first element")
  }
  
  # We use the algorithm detailed in Tibshirani 1994 
  # That uses the Kuhn-Tucker conditions to sequentially find feasible solutions
  
  # First define the starting beta for numerical optimizers
  beta_0<-rep(0,p) 
  
  # Start algorithm by defining delta
  # To do this we need the overall least squares estimate
  # If n>p we can do this explicitly
  # For n<=p we minimize the rss numerically to find a starting value of beta_hat
  
  # Calculate the residual sum of squares and it's derivative for use later
  
  rss<-function(beta) t(y-X%*%beta)%*%(y-X%*%beta)
  rss_deriv<-function(beta) -2*t(X)%*%(y-X%*%beta)
  
  if (n>p){
    # Explicit least-squares solution
    beta_hat<-solve(t(X)%*%X)%*%(t(X) %*% y)
  }else{
    # X'X is not invertable, so use numerical minimization
    beta_hat<-nlm(rss,beta_0)$estimate
  }
  
  # Using this beta_hat we start the algorithm using the K-T conditions
  
  delta=sign(beta_hat)

  # Determine the size of the starting set E
  mod_E=1
  # Calculate the inequality constraint matrix
  G_E<-t(matrix(delta))
  # Inputs for the constrained optimization function have a slightly different form 
  ui<- -G_E
  ci<-rep(-t1,mod_E)
  # Perform a constrained minimization step to attempt to find solution satisfying the lasso constraints
  beta_hat1<-constrOptim(beta_0,rss,rss_deriv,ui,ci,outer.iterations = 200, outer.eps = 1e-10)$par

  # Now start the algorithm to find the constrained solution
  count=1
  while(sum(abs(beta_hat1))>t1){
    # Safety in case the algorthim takes a long time to converge 
    # (guaranteed to converge within 2^p steps, but p could be large)
    if (count>N){
      warning("maximum number of steps reached: consider increasing tol or N")
      break
    }
    # As before calculate delta, then update E and G_E
    delta_new=sign(beta_hat1)
    mod_E<-mod_E+1
    G_E<-rbind(G_E,delta_new)
    # Get constraints in the correct form for constrOptim
    ui<- -G_E
    ci<-rep(-t1,mod_E)

    beta_hat1<-constrOptim(beta_0,rss,rss_deriv,ui,ci)$par
    # The numerical procedures will never shrink coefficients to exactly 0
    # If the absolute value is less than our chosen epsilon we say that is is identically 0
    for (j in 1:p){
      if (abs(beta_hat1[j])<eps){
        beta_hat1[j]=0
      }
    }
    count<-count+1
  }

  for (j in 1:p){
    if (abs(beta_hat1[j])<eps){
      beta_hat1[j]=0
    }
  }
  # Output the value of beta that solves the constrained minimization problem
  beta_hat1
}

