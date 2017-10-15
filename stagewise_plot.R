# for testing stagewise function

set.seed(163)

n <- 14
p <- 9
eps <- 0.05

X <- matrix(rnorm(n*p),nrow=n) #generate some data
y <- X %*% runif(p,-10,10) #generate some y's correlated to X's

mystagewise <- stagewise(X,y,eps,N=850)
B <- mystagewise$coeffs
M <- mystagewise$predicts

par(mfrow=c(2,1),mar=c(4,4,4,2) + 0.1,cex=0.8)

##line graph showing evolution of coefficients beta:
plot(NA, xlim=c(1,ncol(B)+100),ylim=c(-12,12), ylab='betas', xlab='Steps', main='Evolution of betas in stagewise algorithm') 
for (i in 1:p) {
  lines(B[i,],col=i)
}
legend("topright",legend=1:p,col=1:p,lty=1)

##plot of MSE over time:
#MSE <- apply(M,2,function(x){mean((x-y)^2)})
#plot(MSE, type='l',lwd=2, col=2, main='MSE between mu and y over time', ylab='MSE', xlab='Steps')

##plot of AIC over time:
aic <- apply(B,2, function(x){length(x!=0)}) + apply(M,2,function(x){sum((x-y)^2)})
plot(aic, xlim=c(1,ncol(B)+100),ylim=range(aic)+c(-0.5,0.5), type='l',lwd=2, col=2, main='AIC after each step', xlab='steps',ylab="AIC")

X %*% B[,ncol(B)] - y #error after final iteration