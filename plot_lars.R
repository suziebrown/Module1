#' Plot LARS output
#' 
#' Plot method for output from LARS algorithm
#' 
#' 

plot.lars <- function(x, ...) {
  B <- x$beta
  method <- x$method
  
  plot(NA, cex=0.8, xlim=c(1,ncol(B)),ylim=range(B), ylab='betas', xlab='steps', main=sprintf("Evolution of betas in %s", method), ...) 
  for (i in 1:p) {
    lines(B[i,],col=i)
  }
  legend("topleft",legend=1:p,col=1:p,lty=1)
}