#' Cross-validation
#' 
#' @param X matrix of predictor variables
#' @param y vector of response variables
#' @param split_prop proportion to be split into training/test
#' @param n_iterations number random validation procedures
#' @param standardise should data be scaled and centred?
#' @param plot_val Plot the validation error for n_iterations?
#' 
#' 
#' Random Validation works by first of all splitting the whole dataset into trainig and test. We then split the training 
#' again into training2 and validation. We then train on training2 and compute the error on the validation set. Next, we find the 
#' set of betas (number of active betas) which gave us the lowest "validation error", the test set has not been contanimated 
#' unitl now. Lastly we predict on the test set using the betas that gave us minimal "validation error" and compute the "test error". 
#' We repeat n_iterations times.

lars.crossval <- function(X, y, split_prop, n_iterations, standardise, option, plot_val=F){
  # Initialize vectors
  mean_active_set_size = 0
  error_test = c()
  
  for (i in 1:n_iterations){
    # Splitting the data into training, training2 and test sets
    train_index <- sample(1:length(y),length(y)*split_prop)
    test_index <- setdiff(1:length(y), train_index)
    
    train2_index <- sample(train_index, length(train_index)*split_prop)
    val_index <- setdiff(train_index, train2_index)
    
    X_train2 <- X[train2_index,]
    y_train2 <- y[train2_index]
    
    X_val <- X[val_index,]
    y_val <- y[val_index]
    
    X_test <- X[test_index,]
    y_test <- y[test_index]
    
    # Training on the training2 set and computing the validation error
    results <- lars(X_train2, y_train2, option=option, standardise=standardise)
    betas <- results$beta
    predictions_val <- as.matrix(X_val) %*% betas
    
    error_val <- colSums((predictions_val-y_val)**2) + results$t
    
    # Plotting the validation error for n_iterations
    if (plot_val){
      plot(error_val)
    }
    
    # Find the optmial combination of betas
    index_up_to <- which(min(error_val) == error_val)
    mean_active_set_size <- mean_active_set_size + index_up_to/n_iterations
    
    # Computing the test error
    betas_test <- betas[,index_up_to]
    predictions_test <- as.matrix(X_test) %*% betas_test
    error_test <- cbind(error_test, colSums((predictions_test-y_test)**2) + results$t[index_up_to])
  }
  return(list(mean_error_test = mean(error_test), mean_active_set_size = mean_active_set_size))
}
