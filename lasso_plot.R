# For testing lasso function 

n<- 100
p <- 5
X <- matrix(rnorm(n*p),nrow=n) #generate some data
y <- X %*% runif(p,-10,10) #generate some y's correlated to X's

output<-lasso(X,y,1:40,eps=1e-3, standardise=TRUE, intercept=FALSE)
plot(output)



diabetes<-read.csv("diabetes.csv")
dim(diabetes)
head(diabetes)
X<-as.matrix(diabetes[,2:11])
y<-diabetes[,12]

output<-lasso(X,y,seq(from=1,to=4000,by=100),eps=1, standardise=FALSE, intercept=FALSE)
plot(output)


