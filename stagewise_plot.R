# for testing stagewise function

n <- 6
p <- 4
eps <- 0.011

X <- matrix(rnorm(n*p,0,10),nrow=n) #generate some data
y <- X %*% (-2):(p-3) #generate some y's correlated to X's

mystagewise <- stagewise(X,y,eps)
B <- mystagewise$coeffs
M <- mystagewise$predicts

par(mfrow=c(2,1))

##line graph showing evolution of coefficients beta:
plot(NA, xlim=c(1,ncol(B)),ylim=range(B), ylab='betas', xlab='Steps', main='Evolution of betas in stagewise algorithm') 
for (i in 1:p) {
  lines(B[i,],col=i)
}
legend("topright",legend=1:p,col=1:p,lty=1)

##plot of MSE over time:
MSE <- apply(M,2,function(x){mean((x-y)^2)})
plot(MSE, type='l',lwd=2, col=2, main='MSE between mu and y over time', ylab='MSE', xlab='Steps')

X %*% B[,ncol(B)] - y #error after final iteration