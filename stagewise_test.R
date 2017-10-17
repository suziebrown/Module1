## test for stagewise

set.seed(163)

n <- 100
p <- 5
eps <- 0.05

X <- matrix(rnorm(n*p),nrow=n) #generate some data
y <- X %*% runif(p, -2,2) #generate some y's correlated to X's

mystagewise <- stagewise(X,y,eps,N=1000, standardise = F)
plot(mystagewise)

##~~~
## also test Forward:
myforward <- forward(X,y,standardise = F)
plot(myforward)


# 
# B <- mystagewise$beta
# M <- mystagewise$mu
# J <- mystagewise$j
# 
# par(mfrow=c(1,1), cex=0.8)
# 
# plot(NA, xlim=c(1,ncol(B)),ylim=range(B), ylab='betas', xlab='Steps', main='Evolution of betas in stagewise algorithm') 
# for (i in 1:p) {
#   lines(B[i,],col=i)
# }
# legend("topright",legend=1:p,col=1:p,lty=1)

##~~~
## diabetes data

raw <- read.csv("diabetes.csv")
X <- raw[,2:11]
y <- raw[,12]

diab_lasso <- lasso_wrap(X,y,seq(0.05,25,0.05))
plot(diab_lasso)
print(diab_lasso)

diab_stage <- stagewise(X,y,0.1,N=2000)
plot(diab_stage)
print(diab_stage)

diab_forw <- forward(X,y)
plot(diab_forw)
print(diab_forw)
