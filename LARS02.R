#' Least Angle Regression
#' 
#' Implements the LARS algorithm (Efron et. al. (2004) Least angle regression. Annals of Statistics.)
#' 
#' @param X matrix of predictor variables
#' @param y vector of response variables
#' @param standardise should data be scaled and centred?
#' @param t_vec Vector of bounds for the absolute sum of the betas
#' 

lars <-function(X, y, option="lars", t_vec, standardise=T, intercept=F){
  
  n <- nrow(X)
  p <- ncol(X)
  
  ## Edit the data according to options
  if (intercept){
    X <- cbind(rep(1,n), X) #include an intercept term (i.e. a constant covariate)
  }
  if (standardise){
    X <- scale(X)
    y <- y-mean(y)
  }

  ##initialise variables:
  mu <- rep(0,n) #initial estimate is all zeroes
  beta <- rep(0,p) #initial coeficients beta is all zeroes
  M <- mu #matrix containing mu from each step (in columns)
  B <- beta #matrix containing beta from each step (in columns)
  A <- numeric(0) #active set
  J <- NA #vector containing the coordinate changed at each step
  
  ## for first iteration
  corr <- t(X)%*%(y-mu)
  j <- which.max(abs(corr))
  
  ##let's go:
  for (i in 2:(min(n,p))) {
    
    Ac <- setdiff(1:p, A)
    gamma.tmp <- matrix(NA,ncol=length(Ac),nrow=2)
    
    corr <- t(X)%*%(y-mu)
    j <- which.max(abs(corr))
    s_A <- sign(corr[A])
    X_A <- rep(s, each=n) * X[,A]
    G_A <- t(X_A)%*%X_A
    alpha_A <- (sum(G_A))^(-0.5)
    w_A <- alpha_A * rowSums(G_A)
    u_A <- X_A %*% w_A
    a <- t(X)%*%u_A
    
    for (j in 1:length(Ac)) {
      jj <- Ac[j]
      gamma_tmp[,j] <- c((C-corr[jj])/(alpha_A-a[jj]), (C+corr[jj])/(alpha_A+a[jj]))
    }
    j.new <- which.min(gamma_tmp[gamma_tmp > 0])
    
    # ... not finished
  }
  
}