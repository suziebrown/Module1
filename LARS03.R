#' Least Angle Regression
#' 
#' Implements the LARS algorithm (Efron et. al. (2004) Least angle regression. Annals of Statistics.)
#' 
#' @param X matrix of predictor variables
#' @param y vector of response variables
#' @param option which variant of LARS should we run? one of "lars", "lasso", "stagewise"
#' @param standardise should data be scaled and centred?
#' 

lars <-function(X, y, option="lars", standardise=TRUE, intercept=FALSE){

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
  
  ## check inputs are sensible
  if (length(y)!=n){
    stop("number of observations on response not equal to number of observations on predictors")
  }
  if (!(option == "lars" || option == "lasso" || option == "stagewise")){
    stop("invalid option argument, should be one of lars, lasso, stagewise")
  }

  ## Initialise variables
  m <- min(p,n-1) #maximal number of variables in the final active set
  beta <- matrix(0, p,1)
  i <- 0
  mu_old <- matrix(0, n, 1)
  mu <- matrix(0, n, 1)  # Mean Vector
  gamma <- numeric()  # LARS step lengths
  A <- numeric(0) #active set
  Ac <- 1:p #inactive set
  signOK <- 1 #sign check for lasso
  
  ## Let's go
  while (length(A) < m){
    i <- i+1
    
    corr <- t(X)%*%(y-mu) # Current correlation
    C <- max(abs(corr)) # Maximal current absolute correlation
    
    ## For initial step
    if (i == 1){
      j <- which.max(abs(corr))
    } 
    if (signOK){ # Check for LASSO
      A <- c(A,j)
    }
    Ac <- setdiff(1:p, A) # update inactive set 
    
    ## Do all the linear algebra
    s_A <- sign(corr[A])
    X_A <- t(t(X[,A]) * s_A)
    G_A <- t(X_A) %*% X_A
    invG_A <- solve(G_A)
    alpha_A <- 1/sqrt(sum(invG_A)) # normalizing constant
    w_A <- alpha_A[1] * rowSums(invG_A) # Coefficients of equiangular vector u_A
    u_A <- X_A %*% w_A  # Equiangular vector
    a <- t(X) %*% u_A # Angles/correlation between x_j and u_A
    
    # matrix to hold the two possibilities for the minimization to find gamma
    beta_tmp <- matrix(0, p, 1)
    gammaTest <- matrix(0, length(Ac), 2) 
    
    # If we are using all the covariates (last iteration)
    if (length(A) == m){
      gamma[i] <- C/alpha_A   # Move to the least squares projection
    }
    else{
      for (k in 1:length(Ac)){
        jj <- Ac[k]
        gammaTest[k,] <- c((C-corr[jj])/(alpha_A-a[jj]), (C+corr[jj])/(alpha_A+a[jj]))
      }
      gamma[i] <- min(gammaTest[gammaTest>0]) # Take the min over only the positive components
      min_j <- which(gammaTest==min(gammaTest[gammaTest>0]),arr.ind = TRUE)[1] #find ROW of gamma_test holding the min+
      j <- unique(Ac[min_j]) # index to be added to active set
    }
    
    # Update coefficient estimates
    beta_tmp[A] <- t(beta[A,i]) + gamma[i]*w_A*s_A 
    
    # Lasso modification
    if (option=="lasso"){
      signOK <- 1
      # Define the vector d to check for sign changes
      d_A<-s_A*w_A
      gammaTest <- -t(beta[A,i])/d_A
      
      # Look for the first sign change, returning Inf if there are no postive components
      if (all(gammaTest<=0)){
        gamma2 <- Inf 
      }else{
        gamma2 <- min(gammaTest[gammaTest>0]) # Take the min over only the positive components
        min_j <- which(gammaTest==min(gammaTest[gammaTest>0]),arr.ind = TRUE)[2]
      }
      if (gamma2 < gamma[i]){ #The case when sign consistency gets violated
        gamma[i] <- gamma2
        beta_tmp[A] <- beta[A,i] + gamma[i]*d_A    # Correct the coefficients
        beta_tmp[A[unique(min_j)]] <- 0
        A <- A[-unique(min_j)]  # Delete the zero-crossing variable (keep the ordering)
        signOK <- 0
      }
    }
    
    mu <- mu_old + gamma[i] * u_A # Update mean vector
    mu_old <- mu
    beta <- cbind(beta, beta_tmp)  
  }
  
  ## which method name should be returned
  if (option == "lars"){method <- "LARS"}
  if (option == "lasso"){method <- "LARS-Lasso"}
  if (option == "stagewise"){method <- "LARS-Stagewise"}
  
  out <- list(beta = beta, mu = mu, t = colSums(abs(beta)), method=method)
  class(out) <- "lars"
  out
}

