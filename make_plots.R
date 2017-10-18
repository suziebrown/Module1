raw <- read.csv("diabetes.csv")
X <- raw[,2:11]
y <- raw[,12]

par(mfrow=c(2,2))
plot(lars(X,y))
plot(lars(X,y,option='lasso'))
plot(lasso(X,y,seq(1,160,1)))
plot(stagewise(X,y,0.5))
