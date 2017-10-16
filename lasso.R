#' Lasso
#' 
#' Implements the Lasso algorithm
#' 
#' @param X matrix of predictor variables
#' @param y vector of response variables
#' @param t shrinkage parameter
#' 
set.seed(163)
X<-matrix(rnorm(100,mean=1),20,5)
y=X[,1]+2*X[,2]+3*X[,3]+4*X[,4]+5*X[,5]+rnorm(20)

n <- 14
p <- 9
X <- matrix(rnorm(n*p),nrow=n) #generate some data
y <- X %*% runif(p,-10,10) #generate some y's correlated to X's


# Define function to convert binary vectors to decimal

BinToDec <- function(x){
  sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))
}


lasso<-function(X,y,t){
  n <- nrow(X) #number of observations
  p <- ncol(X) #number of covariates
  if (length(y)!=n){stop("number of observations on response not equal to number of observations on predictors")}
  ind <- 1:p
  beta_0<-rep(0,p) # Starting beta for numerical procedures
  
  M <- matrix(NA, nrow=n, ncol=p+1) #matrix containing mu from each step (in columns)
  B <- matrix(NA, nrow=p, ncol=p+1) #matrix containing beta from each step (in columns)
  M[,1] <- rep(0,n) #initial estimate is all zeroes
  B[,1] <- rep(0,p) #initial coeficients beta is all zeroes
  mu <- rep(0,n)

  # We use the algorithm detailed in Tibshirani 1994 
  # That uses the Kuhn-Tucker conditions to sequentially find feasible solutions
  
  # First define delta
  # To do this we need the overall least squares estimate
  # If n>p we can do this explicitly
  # For n<=p we minimize the rss numerically to find a starting value of beta_hat
  
  if (n>p){
    # Explicit least-squares solution
    beta_hat<-solve(t(X)%*%X)%*%(t(X) %*% y)
  }else{
    # X'X is not invertable, so use numerical minimization
    
    rss<-function(beta) t(y-X%*%beta)%*%(y-X%*%beta)
    # Also calculate derivative of rss
    rss_deriv<-function(beta) -2*t(X)%*%(y-X%*%beta)
    beta_hat<-nlm(rss,beta_0)$estimate
  }
  
  # Using this beta_hat we start the algorithm using the K-T conditions
  
  delta=sign(beta_hat)
  
  i_0<-BinToDec(delta)

  E=c(i_0)
  # Need to find the size of E to determine the dimension of the constraint vector
  m<-length(E)
  
  G_E<-t(matrix(delta))
  ui<- -G_E
  ci<-rep(-t,m)
  
  beta_hat1<-constrOptim(beta_0,rss,rss_deriv,ui,ci)$par
  N<-100
  count=1
  while(sum(abs(beta_hat1))>t){
    if (count>N){
      warning("maximum number of steps reached: consider increasing tol or N")
      break
    }
    delta_new=sign(beta_hat1)
    
    i<-BinToDec(delta_new)
    E<-c(E,i)
    G_E<-rbind(G_E,delta_new)
    
    ui<- -G_E
    ci<-rep(-t,m)
    
    beta_hat1<-constrOptim(beta_0,rss,NULL,ui,ci)$par
    for (j in 1:p){
      if (abs(beta_hat1[j])<1e-5){
        beta_hat1[j]=0
      }
    }
    count<-count+1
  }

  for (j in 1:p){
    if (abs(beta_hat1[j])<1e-4){
      beta_hat1[j]=0
    }
  }
  beta_hat1
}

beta_mat<-beta_0
for (t in 1:30){
  beta_new<-lasso(X,y,t)
  beta_mat<-c(beta_mat,beta_new)
}
beta_mat