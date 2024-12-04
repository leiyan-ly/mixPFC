foo1 <- function(i) {
  # split the data
  #set.seed(i)
  #ind_test <- sample(1:N, ntest, replace = FALSE)
  ind_test <- test_ind_mat[i, ] # read saved indices
  x_train <- X[-ind_test, ]
  x_train <- scale(x_train, center = FALSE, scale = TRUE)
  scale_x_train <- attr(x_train, "scaled:scale")
  y_train <- Y[-ind_test]
  x_test <- X[ind_test, ]
  for (ii in seq_len(ncol(x_test))) {
    x_test[, ii] <- x_test[, ii] / scale_x_train[ii] 
  }
  y_test <- Y[ind_test]
  
  # ------ mixPFC ----- #
  mixpfc_mse <- ind_mixpfc <- NA # initialize with NA
  tryCatch({
  # perform estimation using training data
  nfolds <- 5
  foldid <- sample(rep(seq(nfolds), length.out = ntrain))
  # generate intial values
  dcor_vec <- apply(x_train, 2, function(x) dcor(x, y_train))
  select_ind <- tail(order(dcor_vec), d_dcor)
  x_train_cut <- x_train[, select_ind, drop = FALSE]
  pca_fit <- prcomp(x_train_cut, scale. = FALSE, rank. = d_pca)
  x_train_cut_reduced <- pca_fit$x
  # kmeans on x_train_cut_reduced*Y
  xy_train <- x_train_cut_reduced * matrix(y_train, nrow = nrow(x_train_cut_reduced), 
                                           ncol = ncol(x_train_cut_reduced))
  kmeanXY_comp_ind <- kmeans(xy_train, K, nstart = 20)$cluster
  gams_init <- gen_gams_init_K(kmeanXY_comp_ind, true_comp_ind = FALSE)
  # short run of mixPFC on the reduced data
  em_init <- mixPFC_K(x_train_cut_reduced, y_train, FUN = FUN, 
                      n_components = K, gams_init = gams_init,
                      lam_factor = rep(0, K), weight_init = colMeans(gams_init), 
                      center_F_marg = TRUE, compute_loglk = FALSE,
                      max_iter = 15, tol = tol_init)
  gams_init <- em_init$gams 

  mixpfc_fit <- mixPFC_K(x_train, y_train, foldid, FUN = FUN,
                         n_components = K, gams_init = gams_init,
                         lam_factor = lam_factor, weight_init = colMeans(gams_init),
                         center_F_marg = TRUE, compute_loglk = FALSE,
                         print_info = FALSE, max_iter = em_maxiter)
  beta_list <- vector("list", K)
  for (j in seq_len(K)) {
    beta_list[[j]] <- svd(mixpfc_fit$B_mats[[j]])$u[, 1:d_subspace, drop=FALSE]
  }
  
  # estimate the cluster of testing data
  params_list <- list(B_mats = mixpfc_fit$B_mats, 
                      A_mats = mixpfc_fit$A_mats, 
                      Delta_inv = mixpfc_fit$Delta_inv, 
                      weight_pi = mixpfc_fit$weight_pi, 
                      mu = mixpfc_fit$mu)
  gam_test <- est_gam(x_test, y_test, params_list, FUN = FUN)
  
  comp_test_hat <- est_component(gam_test)
  comp_train_hat <- est_component(mixpfc_fit$gams)
  x_test_new <- matrix(NA, nrow = ntest, ncol = d_subspace)
  x_train_new <- matrix(NA, nrow = ntrain, ncol = d_subspace)
  for (j in seq_len(ntest)) {
    k_hat <- comp_test_hat[j]
    x_test_new[j, ] <- x_test[j, ,drop=FALSE] %*% beta_list[[k_hat]]
  }
  
  for (j in seq_len(ntrain)) {
    k_hat <- comp_train_hat[j]
    x_train_new[j, ] <- x_train[j, ,drop=FALSE] %*% beta_list[[k_hat]]
  }
  # compute prediction error for each cluster
  y_test_hat <- rep(NA, ntest)
  for (j in seq_len(K)) {
    train_j <- which(comp_train_hat == j)
    test_j <- which(comp_test_hat == j)
    x_train_new_j <- x_train_new[train_j, , drop=FALSE]
    x_test_new_j <- x_test_new[test_j, , drop=FALSE]
    y_train_j <- y_train[train_j]
    
    lm_fit_j <- lm(y_train_j ~ x_train_new_j)
    coef_lm_j <- matrix(coef(lm_fit_j), ncol = 1)
    y_test_hat[test_j] <- cbind(rep(1, nrow(x_test_new_j)), x_test_new_j) %*% coef_lm_j
    
  }
  
  mixpfc_mse <- mean((y_test - y_test_hat)^2)
  mixpfc_mse_exp <- mean((exp(y_test) - exp(y_test_hat))^2)
  ind_mixpfc <- lapply(mixpfc_fit$B_mats, get_ind_row)
  
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  list(mixpfc_mse = mixpfc_mse, ind_mixpfc = ind_mixpfc,
       mixpfc_mse_exp = mixpfc_mse_exp)
  
}