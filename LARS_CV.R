#' Least Angle Regression
#' 
#' Implements the LARS algorithm (Efron et. al. (2004) Least angle regression. Annals of Statistics.)
#' 
#' @param X matrix of predictor variables
#' @param y vector of response variables
#' @param standardise should data be scaled and centred?
#' @param t_vec Vector of bounds for the absolute sum of the betas
#' 

lars <-function(X, y, option="lars", t_vec, standardise=T){
  
  # option = 'lar' implements the vanilla LAR algorithm (default);
  # option = 'lasso' solves the lasso path with a modified LAR algorithm.
  # t -- a vector of increasing positive real numbers. If given, LARS stops and 
  # returns the solution at t.
  #
  # Output:
  #  A -- a sequence of indices that indicate the order of variable inclusions;
  # beta: history of estimated LARS coefficients;
  # mu -- history of estimated mean vector;
  # C -- history of maximal current absolute corrrelations;
  # c -- history of current corrrelations;
  # gamma: history of LARS step size.
  # Note: history is traced by rows. If t is given, beta is just the
  # estimated coefficient vector at the constraint ||beta||_1 = t.
  #
  # This lars function automatically terminates when the current correlations for inactive set are
  # all zeros. The recovered coefficient vector is the last column of beta 
  # with the *lasso* option. 
  
  
  # Standardising the data 
  if (standardise){
    X = scale(X)
    y = y-mean(y)
  }
  
  eps = 1e-10    # Effective zero
  
  # Still need to work in the Lasso, currently set to LARS
  lasso = F
  
  n <- nrow(X)
  p <- ncol(X)
  
  # Maximal number of variables in the final active set
  m = min(p,n-1)
  
  # Double check what this parameter is doing
  T = length(t)
  
  # Initializing the vectors
  beta = matrix(0, 1, p)
  i = 0
  mu_old = matrix(0, n, 1)
  t_prev = 0
  beta_t = matrix(0, T, p)
  ii = 1
  tt = t_vec
  
  # Mean Vector
  mu = matrix(0, n, 1) 
  
  # LARS step lengths
  gamma = numeric() 
  
  # Active Set
  A = c()
  
  # Not Active set
  Ac = 1:p
  
  # number of variables in the current model
  nVars = 0
  
  ##################################################
  # Double  CHECK  FOR LASSO
  signOK = 1
  ##################################################
  
  # LARS loop
  while (nVars < m){
    i = i+1
    # Current correlation
    corr = t(X)%*%(y-mu) 
    
    # Maximal current absolute correlation
    C = max(abs(corr)) 
    
    # Early stopping criteria (DONT REALLY NEED THIS)
    if (C < eps || length(t)<1){
      break 
    }
    
    # For initial step
    if (i == 1){
      addVar = match(C, abs(corr))
    } 
    
    # Check for LASSO
    if (signOK){
      A = c(A,addVar)
      nVars = nVars+1 # Add one variable to active set
    }
    
    ###############################################################
    # New version
    ###############################################################
    s_A = sign(corr[A])
    #print('THESE ARE THE SA')
    #print(s_A)
    Ones_A = rep(1, length(A))
    
    # Update the Inactive set
    Ac = setdiff(1:p,A)    
    nZeros = length(Ac)
    
    # Gram matrix
    X_A = t(t(X[,A]) * s_A)
    G_A = t(X_A) %*% X_A
    invG_A = solve(G_A)
    
    # Computing the normalizing constant
    A_A = 1/sqrt(sum(invG_A))
    
    # Coefficients of equiangular vector u_A
    w_A = A_A[1] * rowSums(invG_A) 
    
    # Equiangular vector
    u_A = X_A %*% w_A  
    
    # Angles/correlation between x_j and u_A
    a = t(X) %*% u_A 
    
    # matrix to hold the two possibilities for the minimization to find gamma
    beta_tmp = matrix(0, p, 1)
    gammaTest = matrix(0, nZeros, 2) 
    
    # If we are using all the covaraites
    if (nVars == m){
      gamma[i] = C/A_A   # Move to the least squares projection
    }
    else{
      for (j in 1:nZeros){
        jj = Ac[j]
        
        # Computing gamma for the LARS step
        gammaTest[j,] = c((C-corr[jj])/(A_A-a[jj]), (C+corr[jj])/(A_A+a[jj]))
      }
      
      # Take the min over only the positive components
      gamma[i] = min(gammaTest[gammaTest>0]) 
      min_j = which(gammaTest==min(gammaTest[gammaTest>0]),arr.ind = TRUE)[1]
      
      # index to add into the Active set
      addVar = unique(Ac[min_j])
    }
    
    # Update coefficient estimates
    beta_tmp[A] = beta[i,A] + gamma[i]*w_A*s_A 
    
    # Check the sign feasibility of lasso
    if (lasso){
      signOK = 1
      gammaTest = -t(beta[i,A])/w_A
      gamma2 = min(gammaTest[gammaTest>0]) # Take the min over only the positive components
      min_j = which(gammaTest==min(gammaTest[gammaTest>0]),arr.ind = TRUE)[1]
      if (gamma2 < gamma[i]){ #The case when sign consistency gets violated
        gamma[i] = gamma2
        beta_tmp[A] = t(beta[i,A]) + gamma[i]*w_A    # Correct the coefficients
        beta_tmp[A[unique(min_j)]] = 0
        A[unique(min_j)] = numeric()  # Delete the zero-crossing variable (keep the ordering)
        nVars = nVars-1
        signOK = 0
      }
    }
    
    
    # # Need to reread this works for t as 1-D vector not sure how it is supposed to work in more than 1
    # #print(t_vec[1])
    # if (Inf != t_vec[1]){
    #   t_now = sum(abs(beta_tmp[A]))
    #   if (t_prev < t_vec[1] && t_now >= t_vec[1]){
    #     beta_t[ii,A] = beta[i,A] + A_A %*% (t_vec[1]-t_prev) %*% t(w_A)    # Compute coefficient estimates corresponding to a specific t
    #     #t_vec = t_vec[-1]
    #     ii = ii+1
    #   }
    #   t_prev = t_now
    # }
    
    mu = mu_old + gamma[i] * u_A # Update mean vector
    mu_old = mu
    beta = rbind(beta, t(beta_tmp))  
  }
  
  list(beta = t(beta), J = A, mu = mu, C = C, c = c, gamma = gamma, t = colSums(abs(t(beta))), method = "LARS")
}


Data <- read.csv(file="diabetes.csv", header=TRUE, sep=",")
Data <- Data[,-1]
X <- Data[, 1:10]
y <- Data[, 11]
e <- rep(0, 11)
#X <- scale(X)
y <- y-mean(y)


 # library(MASS)
 # data('Boston')
 # X <- as.matrix(Boston[,1:13])
 # y <- Boston[,14]
 # 
 # X <- scale(X)
 # y <- y-mean(y)



Cross_Validation_LARS <- function(X, y, split_prop, n_iterations, standardise, option){
  mean_active_set_size = 0
  error_test = c()
  for (i in 1:n_iterations){
    train_index <- sample(1:length(y),length(y)*split_prop)
    test_index <- setdiff(1:length(y), train_index)
    
    train2_index <- sample(train_index, length(train_index)*split_prop)
    val_index <- setdiff(train_index, train2_index)
    
    X_train2 <- X[train2_index,]
    y_train2 <- y[train2_index]
    
    X_val <- X[val_index,]
    y_val <- y[val_index]
    
    X_test <- X[test_index,]
    y_test <- y[test_index]
    
    results <- lars(X_train2, y_train2, option="lars", t_vec= c(10) , standardise=standardise)
    betas <- results$beta
    predictions_val <- as.matrix(X_val) %*% betas
    
    error_val <- colSums((predictions_val-y_val)**2) + results$t
    
    plot(error_val)
    # e <- e + error_val/100 
    
    index_up_to <- which(min(error_val) == error_val)
    mean_active_set_size <- mean_active_set_size + index_up_to/n_iterations
    Active_set <- results$J
    
    betas_test <- betas[,index_up_to]
    predictions_test <- as.matrix(X_test) %*% betas_test
    # Inactive <- setdiff(1:length(results$t), Active_set[1:index_up_to])
    error_test <- cbind(error_test, colSums((predictions_test-y_test)**2) + results$t[index_up_to])
  }
  
  return(list(error_test = as.vector(error_test), mean_active_set_size = mean_active_set_size) )
}


res <- Cross_Validation_LARS(X = X, y=y , split_prop = 0.9, n_iterations = 100, standardise = T, option = 'lars')

res$error_test

# class(results) <- "lars"
# plot(results)
# mean(error_test/length(y_test))
# 
# plot(as.vector(error_test))
# error_test
# results$J
# mean_active_set_size
