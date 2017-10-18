# Test data that should cause problems for lasso/lars

# 1) Use data taht violates the irresentable condition
# ie that significant variables cannot be correlated with insignificant ones

library(mvtnorm)
library(corrplot)
library(glmnet)
library(clusterGeneration)

k=10 # = Number of Candidate Variables
p=5 # = Number of Relevant Variables
N=500 # = Number of observations
betas=(-1)^(1:p) # = Values for beta
set.seed(12345) # = Seed for replication
sigma1=genPositiveDefMat(k,"unifcorrmat")$Sigma # = Sigma1 violates the irc
sigma2=sigma1 # = Sigma2 satisfies the irc
sigma2[(p+1):k,1:p]=0
sigma2[1:p,(p+1):k]=0

# = Verify the irrepresentable condition
irc1=sort(abs(sigma1[(p+1):k,1:p]%*%solve(sigma1[1:p,1:p])%*%sign(betas)))
irc2=sort(abs(sigma2[(p+1):k,1:p]%*%solve(sigma2[1:p,1:p])%*%sign(betas)))
c(max(irc1),max(irc2))


# = Have a look at the correlation matrices
par(mfrow=c(1,2))
corrplot(cov2cor(sigma1))
corrplot(cov2cor(sigma2))

X1=rmvnorm(N,sigma = sigma1) # = Variables for the design that violates the IRC = #
X2=rmvnorm(N,sigma = sigma2) # = Variables for the design that satisfies the IRC = #

#X1<- scale(X1) #centre and normalise X1
#X2 <- scale(X2) #centre and normalise X2

e=rnorm(N) # = Error = #
y1=X1[,1:p]%*%betas+e # = Generate y for design 1 = #
y2=X2[,1:p]%*%betas+e # = Generate y for design 2 = #

#y1 <- y1 - mean(y1) #centre y
#y2 <- y2 - mean(y1) #centre y

lasso1=glmnet(X1,y1,nlambda = 100) # = Estimation for design 1 = #
lasso2=glmnet(X2,y2,nlambda = 100) # = Estimation for design 2 = #

## == Regularization path == ##
par(mfrow=c(1,2))
l1=log(lasso1$lambda)
matplot(as.matrix(l1),t(coef(lasso1)[-1,]),type="l",lty=1,col=c(rep(1,9),2),ylab="coef",xlab="log(lambda)",main="Violates IRC")
l2=log(lasso2$lambda)
matplot(as.matrix(l2),t(coef(lasso2)[-1,]),type="l",lty=1,col=c(rep(1,9),2),ylab="coef",xlab="log(lambda)",main="Satisfies IRC")

lasso_lars_out1 <- lars(X = X1, y = y1, option='lasso',  t_vec = c(1,2,3,4,5,6,7,8), standardise = F)
lasso_lars_out2 <- lars(X = X2, y = y2, option='lasso',  t_vec = c(1,2,3,4,5,6,7,8), standardise = F)

class(lasso_lars_out1) <- "lars"
class(lasso_lars_out2) <- "lars"

plot(lasso_lars_out1)
plot(lasso_lars_out2)

lasso_out1<-lasso(X1,y1,1:10,eps=1e-2, standardise=FALSE, intercept=FALSE)
lasso_out2<-lasso(X2,y2,1:10,eps=1e-2, standardise=FALSE, intercept=FALSE)

plot(lasso_out1)
plot(lasso_out2)



lars_out1 <- lars(X = X1, y = y1, t_vec = c(1,2,3,4,5), standardise = F)
lars_out2 <- lars(X = X2, y = y2, t_vec = c(1,2,3,4,5), standardise = F)

class(lars_out1) <- "lars"
class(lars_out2) <- "lars"

# plot(lars_out1)
# plot(lars_out2)
pdf("break_plots.pdf",width=8, height=5)
par(mfrow=c(1,2))
op <- par(cex = 1)

x<- lars_out1
B <- x$beta
method <- x$method
t <- x$t
p <- nrow(B)

plot(NA, cex=0.8, xlim=range(t),ylim=range(B), ylab='betas', xlab='t', main="Evolution of betas using Sigma_1") 
for (i in 1:5) {
  lines(t, B[i,],col=i)
}
for (i in 2:p) {
  lines(t, B[i,],col=i,lty=2)
}
op <- par(cex = 0.5)

legend("topleft",legend=c("X_1","X_2","X_3","X_4","X_5","X_6","X_7","X_8","X_9","X_10"),col=1:p,lty=c(rep(1,5),rep(2,p-5)))

op <- par(cex = 1)

x<- lars_out2
B <- x$beta
method <- x$method
t <- x$t
p <- nrow(B)

plot(NA, cex=0.8, xlim=range(t),ylim=range(B), ylab='betas', xlab='t', main="Evolution of betas using Sigma_2") 
for (i in 1:5) {
  lines(t, B[i,],col=i)
}
for (i in 2:p) {
  lines(t, B[i,],col=i,lty=2)
}
op <- par(cex = 0.5)

legend("topleft",legend=c("X_1","X_2","X_3","X_4","X_5","X_6","X_7","X_8","X_9","X_10"),col=1:p,lty=c(rep(1,5),rep(2,p-5)))

dev.off()


