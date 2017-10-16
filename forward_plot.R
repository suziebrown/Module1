# plots from forward function

set.seed(163)

n <- 14
p <- 9
X <- matrix(rnorm(n*p),nrow=n) #generate some data
y <- X %*% runif(p,-10,10) #generate some y's correlated to X's

myforward <- forward(X,y)
B <- myforward$coeffs
M <- myforward$predicts

par(mfrow=c(2,1),mar=c(4,4,4,2) + 0.1,cex=0.8)

##line graph showing evolution of coefficients beta:
plot(NA, xlim=c(1,ncol(B)+1),ylim=c(-12,12), ylab='betas', xlab='steps', main='Evolution of betas in forward algorithm') 
for (i in 1:p) {
  points(B[i,],col=i, type='b', pch=16)
}
legend("topright",legend=1:p,col=1:p,lty=1, pch=16)

##plot of AIC/MSE over time:
aic <- 2*(0:(ncol(B)-1)) + apply(M,2,function(x){n*log(sum((x-y)^2))})
plot(aic, xlim=c(1,ncol(B)+1),ylim=range(aic)+c(-0.5,0.5), type='l',lwd=2, col=2, main='AIC after each step', xlab='steps',ylab="AIC")
#MSE <- apply(M,2,function(x){mean((x-y)^2)})
#lines(MSE, type='l',lwd=2, col=4)
#legend("topright",c("MSE","AIC"), col=c(4,2),lwd=2)
