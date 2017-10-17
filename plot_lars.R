#' Plot LARS output
#' 
#' Plot method for output from LARS algorithm
#' 
#' @param x an object of class "lars"
#' 

plot.lars <- function(x, ...) {
  B <- x$beta
  method <- x$method
  t <- x$t
  p <- nrow(B)
  
  plot(NA, cex=0.8, xlim=range(t),ylim=range(B), ylab='betas', xlab='t', main=sprintf("Evolution of betas in %s", method), ...) 
  for (i in 1:p) {
    lines(t, B[i,],col=i)
  }
  legend("topleft",legend=1:p,col=1:p,lty=1)
}

