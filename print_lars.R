#' Summarise LARS output
#' 
#' Summary method for output from LARS algorithm
#' 
#' @param x an object of class "lars"
#' 

print.lars <- function(x, ...) {
  B <- x$beta
  method <- x$method
  t <- x$t
  J <- x$j
  
  cat("Linear regression model selection using ", method, "\n")
  
  if (method=="LARS" || method=="Forward") {
    cat(ncol(B), " models proposed \n")
    cat("\n coefficients for each model: \n")
    print(B)
    cat("\n L1 penalty for each model: \n")
    print(t)
  }
  
  if (method=="Stagewise"){
    
    A_old <- numeric(0)
    use <- numeric(0)
    for (i in 1:length(J)) {
      A_new <- unique(J[1:i])
      if (length(setdiff(A_new, A_old))!=0) {
        use <- c(use, i-1)
      }
      A_old <- A_new
    }
    use <- c(use, length(J))
    
    cat(length(use)-1, " models proposed \n")
    cat("\n coefficients for each model: \n")
    print(B[,use])
    cat("\n L1 penalty for each model: \n")
    print(t[use])
  }
  
}