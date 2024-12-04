# the function to compute predictive error on the testing set
pred_err <- function(x_train, y_train, x_test, y_test) {
  # x_train, y_train: training data
  # x_test, y_test: testing data
  # both x_train and x_test should be matrix, even if they only have one column
  # ----- linear regression ----- #
  lm_fit <- lm(y_train ~ x_train)
  coef_lm <- matrix(coef(lm_fit), ncol = 1) # remember it includes the intercept
  y_test_hat <- cbind(rep(1, nrow(x_test)), x_test) %*% coef_lm
  
  log_mse <- mean((y_test - y_test_hat)^2)
  mse <- mean((exp(y_test) - exp(y_test_hat))^2)

  c(log_mse, mse)

}

compute_mse <- function(y_test, y_test_hat) {
  log_mse <- mean((y_test - y_test_hat)^2)
  mse <- mean((exp(y_test) - exp(y_test_hat))^2)
  
  c(log_mse, mse)
  
}

foo2 <- function(i) {
  ind_test <- test_ind_mat[i, ] # read saved indices
  x_train <- X[-ind_test, ]
  x_train <- scale(x_train, center = FALSE, scale = TRUE)
  scale_x_train <- attr(x_train, "scaled:scale")
  y_train <- Y[-ind_test]
  x_test <- X[ind_test, ]
  for (ii in seq_len(ncol(x_test))) {
    #cat(ii, "\n")
    x_test[, ii] <- x_test[, ii] / scale_x_train[ii] 
  }
  y_test <- Y[ind_test]
  # calculate the prediction error on the testing data
  nfolds <- 5
  foldid <- sample(rep(seq(nfolds), length.out = ntrain))
  ### treat as one cluster ###
  # ----- seasSIR ----- #
  seassir_fit <- cv.seas(x_train, y_train, H = H, foldid = foldid, d = d_subspace,
                         lam1_fac=seq(1.5,0.8, length.out = 10), lam2 = 0, 
                         type = 'sir') 
  beta_seassir <- seassir_fit$beta
  if (is.null(beta_seassir) || (seassir_fit$rank == 0)) {
    seassir_mse <- NA
    ind_seassir <- NA
  } else {
    x_train_new <- x_train %*% beta_seassir                              # Reduce the training set.
    x_test_new <- x_test %*% beta_seassir                                # Reduce the test set.
    seassir_mse <- pred_err(x_train_new, y_train, x_test_new, y_test)
    ind_seassir <- get_ind_row(beta_seassir)
  }
  
  # ----- seasPFC ----- #
  seaspfc_fit <- cv.seas(x_train, y_train, foldid = foldid, d = d_subspace,
                         lam1_fac=seq(1.5,0.8, length.out = 10), lam2 = 0, 
                         type = 'pfc', FUN = FUN) 
  beta_seaspfc <- seaspfc_fit$beta
  if (is.null(beta_seaspfc) || (seaspfc_fit$rank == 0)) {
    seaspfc_mse <- NA
    ind_seaspfc <- NA
  } else {
    x_train_new <- x_train %*% beta_seaspfc                              # Reduce the training set.
    x_test_new <- x_test %*% beta_seaspfc                                # Reduce the test set.
    seaspfc_mse <- pred_err(x_train_new, y_train, x_test_new, y_test)
    ind_seaspfc <- get_ind_row(beta_seaspfc)
  }
  # ----- lassoSIR ----- #
  lassosir_fit <- LassoSIR_revised(x_train, y_train, H = H, 
                                   foldid = foldid, no.dim = d_subspace) 
  beta_lassosir <- lassosir_fit$beta
  if (is.null(beta_lassosir) || (lassosir_fit$no.dim == 0)) {
    lassosir_mse <- NA
    ind_lassosir <- NA
  } else {
    x_train_new <- x_train %*% beta_lassosir
    x_test_new <- x_test %*% beta_lassosir
    lassosir_mse <- pred_err(x_train_new, y_train, x_test_new, y_test)
    ind_lassosir <- get_ind_row(beta_lassosir)
  }
  
  # ----- lasso ----- #
  cv_lasso <- cv.glmnet(x_train, y_train, nfolds = 5)
  y_test_hat <- predict(cv_lasso, newx = x_test, s = "lambda.min")
  
  lasso_beta <- coef(cv_lasso, s = "lambda.min")[-1]
  ind_lasso <- which(lasso_beta != 0) 
  lasso_mse <- compute_mse(y_test, y_test_hat)

  # ----- lasso no intercept ----- #
  x_train_1 <- cbind(rep(1, nrow(x_train)), x_train)
  cv_lasso_no_intercept <- cv.glmnet(x_train_1, y_train, nfolds = 5, intercept = FALSE)
  x_test_1 <- cbind(rep(1, nrow(x_test)), x_test)
  y_test_hat_1 <- predict(cv_lasso_no_intercept, newx = x_test_1, s = "lambda.min")
  
  lasso_beta_1 <- coef(cv_lasso_no_intercept, s = "lambda.min")[-1]
  ind_lasso_1 <- which(lasso_beta_1 != 0) 
  lasso_mse_1 <- compute_mse(y_test, y_test_hat_1)


  # ----- refit lasso ----- #
  lasso_refit <- lm(y_train ~ x_train[, ind_lasso, drop=FALSE])
  coef_lassorefit <- matrix(coef(lasso_refit), ncol = 1)
  y_test_hat <- cbind(rep(1, nrow(x_test)), x_test[, ind_lasso, drop=FALSE]) %*% coef_lassorefit
  lasso_refit_mse <- compute_mse(y_test, y_test_hat)
  
  ### collect results 
  pred_mse <- cbind(seassir_mse, seaspfc_mse, lassosir_mse, lasso_mse, lasso_refit_mse, lasso_mse_1)
  colnames(pred_mse) <- c("seasSIR", "seasPFC", "lassoSIR", "Lasso", "LassoRefit", "LassoOne")
  sele_var <- list(ind_seassir, ind_seaspfc, ind_lassosir, ind_lasso, ind_lasso, ind_lasso_1)
  names(sele_var) <- c("seasSIR", "seasPFC", "lassoSIR", "Lasso", "LassoRefit", "LassoOne")
  
  list(pred_mse = pred_mse, sele_var = sele_var)
}

