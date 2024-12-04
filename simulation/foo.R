foo <- function(i, em_maxiter = 60) {
  set.seed(i)
  est_comp_hat <- vector("list", n_Delta)
  for (j in j_seq) {
    #cat("simulation setting: ", K, 
    #    ", Delta&Gamma: ", j, ", rep: ", i, "\n")
    
    Delta <- Delta_list[[j]]
    Gamma <- Gamma_list[[j]]
    
    
    #train <- data_gen_K_random(N, p, Delta, Gamma, weight, 
    #                           center_F_marg = TRUE, K = K)
    
    ### read data
    file_path <- paste0(save_data_path, cov_path[j], '/', 'train_', i, '.mat')
    train <- readMat(file_path)[[1]]
    
    x_train <- train[[1]]
    y_train <- train[[2]]
    component_ind <- train[[3]]
    ###
    
    nfolds <- 5
    foldid <- sample(rep(seq(nfolds), length.out = N))
    
    ### generate intial values
    dcor_vec <- apply(x_train, 2, function(x) dcor(x, y_train))
    select_ind <- tail(order(dcor_vec), d_dcor[j])
    x_train_cut <- x_train[, select_ind, drop = FALSE]
    pca_fit <- prcomp(x_train_cut, scale. = FALSE, rank. = d_pca[j])
    x_train_cut_reduced <- pca_fit$x
    # kmeans on x_train_cut_reduced*Y
    xy_train <- x_train_cut_reduced * matrix(y_train, nrow = nrow(x_train_cut_reduced), ncol = ncol(x_train_cut_reduced))
    kmeanXY_comp_ind <- kmeans(xy_train, K, nstart = 20)$cluster
    er_kmeanX_reduced_Y[j] <- error_rate_K(kmeanXY_comp_ind, component_ind, K)$error_rate
    gams_init <- gen_gams_init_K(kmeanXY_comp_ind, true_comp_ind = FALSE)
    # short run of mixPFC on the reduced data
    em_init <- mixPFC_K(x_train_cut_reduced, y_train, FUN = FUN, 
                        n_components = K, gams_init = gams_init,
                        lam_factor = rep(0, K), weight_init = colMeans(gams_init), 
                        center_F_marg = TRUE, compute_loglk = FALSE,
                        max_iter = 15, tol = tol_init[j])
    gams_init <- em_init$gams
    er_init[j] <- error_rate_K(est_component(gams_init), component_ind, K)$error_rate
    

    # run mixpfc with gams_init
    mixpfc_fit <- mixPFC_K(x_train, y_train, foldid, FUN = FUN,
                           n_components = K, gams_init = gams_init,
                           lam_factor = lam_factor[[j]], weight_init = colMeans(gams_init),
                           center_F_marg = TRUE, compute_loglk = FALSE,
                           print_info = FALSE, max_iter = em_maxiter)
    
    # eval results
    component_ind_hat <- est_component(mixpfc_fit$gams)
    
    # do one more e-step
    post_prob <- compute_post_K(x_train, y_train, mixpfc_fit$B_mats,
                                mixpfc_fit$weight_pi, mixpfc_fit$gams, FUN = FUN,
                                d_vec = d_vec, center_F_marg = TRUE)
    component_ind_hat2 <- est_component(post_prob)
    est_comp_hat[[j]] <- component_ind_hat2
    
    er_and_ord <- error_rate_K(component_ind_hat, component_ind, K) # error rate and the best permutation of cluster
    # store the error rates
    er_rate[j] <- er_and_ord$error_rate
    er_rate2[j] <- error_rate_K(component_ind_hat2, component_ind, K)$error_rate
    # error_rate_each_cls(component_ind_hat, component_ind)
    
    B_hat <- mixpfc_fit$B_mats
    beta_hat <- vector("list", length = K)
    for (k in seq_len(K)) {
      beta_hat[[k]] <- svd(B_hat[[k]])$u[, 1:d_vec[k], drop = FALSE]
    }
    
    # store the subspace distance
    subspace_dist_and_ord <- sub_dist_error_K(beta, beta_hat) # subspace distance and the best permutation of beta_hat
    subspace_dist[j, ] <- subspace_dist_and_ord$dist_error
    
    if (!all(er_and_ord$pred_K == subspace_dist_and_ord$matched_perm)) {
      cat("simulation: ", i, "Delta: ", j, "\n")
      cat("the permutation of error rate: ", er_and_ord$pred_K, 
          "| the permutation of subspace distance: ", subspace_dist_and_ord$matched_perm, "\n")
      cat("error rate: ", er_rate[j], 
          "| subspace distance: ", subspace_dist[j, ], "\n")
    }
    # we'll stick with the permutation given by subspace distance
    row_ind_hat <- lapply(B_hat[subspace_dist_and_ord$matched_perm], get_ind_row)
    for (kk in seq_len(K)) {
      TPR[j, kk] <- compute_TPR(row_ind_hat[[kk]], ind_row[[kk]])
      FPR[j, kk] <- compute_FPR(row_ind_hat[[kk]], ind_row[[kk]], p)
    }
    row_ind_hat_union <- Reduce(union, row_ind_hat)
    TPR_union[j] <- compute_TPR(row_ind_hat_union, ind_row_union)
    FPR_union[j] <- compute_FPR(row_ind_hat_union, ind_row_union, p)
    # store estimated weights
    pi_hat[j, ] <- colMeans(mixpfc_fit$gams)
    pi_hat2[j, ] <- colMeans(post_prob)
    
    #cat("TPR of Delta ", j, ": ", TPR[j, ], "\n")
    #cat("FPR of Delta ", j, ": ", FPR[j, ], "\n")
    
    ##### compute oracle  ######
    for (k in seq_len(K)) {
      x_k <- x_train[component_ind == k, ind_row[[k]]]
      y_k <- y_train[component_ind == k]
      MU_k <- MU(x_k, y_k, type = "pfc", FUN = FUN)
      fit_k <- solve(MU_k$M) %*% MU_k$U
      subspace_dist_ora[j, k] <- subspace(beta[[k]][ind_row[[k]], , drop=FALSE], svd(fit_k)$u[, 1:d_vec[k], drop=FALSE])
    }
    # compute the error rate
    if (all(d_vec == 1)) {
      etas <- vector("list", K)
      for (k in seq_len(K)) {
        etas[[k]] <- matrix(c(1, 0.3), nrow = 1)
      }
      FUN_true <- function(y) c(y, abs(y))
      
    } else if (all(d_vec == 2)) {
      etas <- vector("list", K)
      for (k in seq_len(K)) {
        etas[[k]] <- diag(2)
      }
      FUN_true <- function(y) c(y, abs(y))
    }
    
    post_prob_ora <- compute_post_ora_K(x_train, y_train, weight, Gamma, Delta, etas, FUN = FUN_true)
    comp_hat_ora <- est_component(post_prob_ora)
    er_rate_ora[j] <- sum(comp_hat_ora != component_ind) / N
    
  }
  list(er_rate = er_rate, er_rate2 = er_rate2, subspace_dist = subspace_dist,
       pi_hat = pi_hat, pi_hat2 = pi_hat2, subspace_dist_ora = subspace_dist_ora,
       er_rate_ora = er_rate_ora, TPR = TPR, FPR = FPR, TPR_union = TPR_union, 
       FPR_union = FPR_union, er_init = er_init, er_kmeanX_reduced_Y = er_kmeanX_reduced_Y,
       est_comp_hat = est_comp_hat)
}

# save the estimated labels obtained using mixPFC
save_comp_hat <- function(sim_result, path) {
  # === Input
  # sim_result: a list with length n_rep, number of repitions
  # each element is a list, consisting of estimation results for each dataset
  # path: file path to save the estimated labels
  # === Output:
  # this function saves the estimated labels to a given path, no output
  # the saved file contains a list variable, i-th element is the estimated
  # labels for i-th dataset
  n_sim <- length(sim_result)
  comp_hat <- vector("list", n_sim)
  for (i in seq_len(n_sim)) { 
    comp_hat[[i]] <- sim_result[[i]]$est_comp_hat
  }
  
  # save comp_hat
  save(comp_hat, file = path)
}

