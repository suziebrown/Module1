# Testing the algorithm using the Diabetes dataset
set.seed(100)
Data <- read.csv(file="diabetes.csv", header=TRUE, sep=",")
Data <- Data[,-1]
X <- Data[, 1:10]
y <- Data[, 11]
X <- scale(X)
y <- y-mean(y)

res <- lars.crossval(X, y, 0.9, 10, T, 'lasso')
res$mean_active_set_size
res$mean_error_test