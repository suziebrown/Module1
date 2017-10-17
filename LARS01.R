#' Least Angle Regression
#' 
#' Implements the LARS algorithm (Efron et. al. (2004) Least angle regression. Annals of Statistics.)
#' 
#' @param X matrix of predictor variables
#' @param Y vector of response variables
#' @param standardize should data be scaled and centred?
#' @param t_vec Vector of bounds for the absolute sum of the betas

lars <-function(X, Y, option, t_vec, standardize=T){

# Least Angle Regression (LAR) algorithm.
# Ref: Efron et. al. (2004) Least angle regression. Annals of Statistics.
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
# Remarks:
  # 1. LARS is originally proposed to estimate a sparse coefficient vector in
# a noisy over-determined linear system. LARS outputs estimates for all
# shrinkage/constraint parameters (homotopy).
#
# 2. LARS is well suited for Basis Pursuit (BP) purpose in the real case. This lars function
# automatically terminates when the current correlations for inactive set are
# all zeros. The recovered coefficient vector is the last column of beta 
# with the *lasso* option. Hence, this function provides a fast and 
# efficient solution for the ell_1 minimization problem. 
# Ref: Donoho and Tsaig (2006). Fast solution of ell_1 norm minimization problems when the solution may be sparse.

# if nargin < 5, standardize = true; end
# if nargin < 4, t = Inf; end
# if nargin < 3, option = 'lar'; end

# if strcmpi(option, 'lasso'), lasso = 1; else, lasso = 0; end

eps = 1e-10    # Effective zero

# Still need to work in the Lasso, currently set to LARS
lasso = F

n <- nrow(X)
p <- ncol(X)

# Standardizing the data 
if (standardize){
  X = scale(X)
  Y = Y-mean(Y)
}

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
  corr = t(X)%*%(Y-mu) 
  
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
  Ones_A = rep(1, length(A))
  
  # Update the Inactive set
  Ac = setdiff(1:p,A)    
  nZeros = length(Ac)
  
  # Gram matrix
  X_A = t(t(X[,A]) * s_A)
  G_A = t(X_A) %*% X_A
  invG_A = solve(G_A)
  
  # Computing the normalizing constant
  A_A = 1/sqrt(t(Ones_A) %*% invG_A %*% Ones_A)
  
  # Coefficients of equiangular vector u_A
  w_A = A_A[1] * invG_A %*% Ones_A 
  
  # Equiangular vector
  u_A = X_A %*% w_A  
  
  # Angles/correlation between x_j and u_A
  a = t(X) %*% u_A 
  
  # matrix to hold the two possibilities for the minimization to find gamma
  beta_tmp = matrix(0, p, 1)
  gammaTest = matrix(0, nZeros, 2) 
  
  ###############################################################
  
  # LOOK AT THIS AGAIN
  # This method of calculating u_A is different to the paper, 
  # maybe should revert back to paper method - Cian
  #s_A = sign(corr[A])
  #Ac = setdiff(1:p,A)    # Inactive set
  #nZeros = length(Ac)
  #X_A = X[,A] 
  #G_A = t(X_A)*X_A # Gram matrix
  #invG_A = solve(G_A)
  #L_A = 1/sqrt(t(s_A)%*%invG_A%*%s_A)
  #w_A = L_A%*%invG_A*s_A   # Coefficients of equiangular vector u_A
  #u_A = X_A %*% w_A  # Equiangular vector
  #a = t(X) %*% u_A # Angles between x_j and u_A
  #beta_tmp = matrix(0, p, 1)
  #gammaTest = matrix(0, nZeros, 2) # matrix to hold the two possibilities for the minimization to find gamma
  ###################################################################################################
  ###################################################################################################
  
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
  beta_tmp[A] = beta[i,A] + gamma[i]*w_A   
  
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
  
  
  # Need to reread this works for t as 1-D vector not sure how it is supposed to work in more than 1
  print(t_vec[1])
  if (Inf != t_vec[1]){
      t_now = sqrt(sum(beta_tmp[A]**2))
      if (t_prev < t_vec[1] && t_now >= t_vec[1]){
        beta_t[ii,A] = beta[i,A] + A_A %*% (t_vec[1]-t_prev) %*% t(w_A)    # Compute coefficient estimates corresponding to a specific t
        #t_vec = t_vec[-1]
        ii = ii+1
      }
   t_prev = t_now
  }
  mu = mu_old + gamma[i] * u_A # Update mean vector
  mu_old = mu
  beta = rbind(beta, t(beta_tmp))  
}

list(beta = beta, A = A, mu = mu, C = C, c = c, gamma = gamma)
}



Data <- read.csv(file="diabetes.csv", header=TRUE, sep=",")
Data <- Data[,-1]
X <- Data[, 1:10]
Y <- Data[, 11]


t_max = 1000
t_vector <- 1:t_max
results_beta <- vector("list", t_max)
#results_beta = rep(0,10)
for (t in t_vector){
  results <- lars(X = X, Y = Y, t_vec = c(t), standardize = T)
  results_beta[[t]] <- results$beta
}
results <- lars(X = X, Y = Y, t_vec = c(10,20,30,40,50,60), standardize = T, method = "LARS")


betas <- t(results)
class(betas) <- "lars"
plot(betas)

##############################################################################
# MAYBE GO OVER THIS TOGETHER CAN T FIGURE OUT WHAT IS BEING DONE HERE
##############################################################################
# if (1 < ii){
#   noCons = (tt > norm(beta_tmp,1))
#   if (0 < sum(noCons)){
#     ########################################################################
#     # Check repmat function
#     ########################################################################
#     beta_t(noCons, ) = repmat(t(beta_tmp),sum(noCons),1)
#   }
#   beta = beta_t
# }

##############################################################################
##############################################################################

##############################################################################
# MOST PROBABLY NOT NEEDED 
##############################################################################

# Normalize columns of X to have mean zero and length one.
#sX_f <-  function (X){
#  return(scale(X))
#}

#n = dim(X)[1]
#p = dim(X)[2]

#sX = X-repmat(mean(X),n,1)
#sX = sX*diag(1./sqrt(ones(1,n)*sX.^2))

#sX <- sX_f(X)

############################################################
# Check the Minplus function
############################################################

# Find the minimum and its index over the (strictly) positive part of X
# matrix
#function [m, I, J] = minplus(X)

# Remove complex elements and reset to Inf
#[I,J] = find(0~=imag(X))
#for i = 1:length(I),
#  X(I(i),J(i)) = Inf
#end

#X(X<=0) = Inf
#m = min(min(X))
#[I,J] = find(X==m)