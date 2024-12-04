# this file implements mixture PFC with isotonic covariance matrix and eta

# E-step: Compute the expected value of the latent variables
e_step_iso_K <- function(X, Fmat, params, center_F_marg = FALSE, compute_loglk = TRUE) {
  # Implement the E-step to compute expected values of latent variables
  # Inputs:
  # ======
  # X: n*p design matrix
  # Fmat: F matrix of dimension n*q
  # params: a list of params, consists of 
  #        Gam_mats, weight_pi, gams
  #        sigma2, eta_mats
  # center_F_marg: boolean, if to center F(y) marginally
  # compute_loglk: boolean, if to compute loglikelihood
  # Outputs:
  # =======
  # gams: updated responsibilities
  # loglk: log likelihood of current estimates
  
  Gam_mats <- params$Gam_mats
  eta_mats <- params$eta_mats
  gams <- params$gams
  weight_pi <- params$weight_pi
  sigma2 <- params$sigma2
  

  n <- nrow(gams) # number of samples
  p <- ncol(X) # number of variables
  n_components <- ncol(gams) # number of clusters
  
  X_c_mats <- Fmat_c_mats <- vector("list", length = n_components)
  # cut the X_c matrix to only keep the selected variables
  # X_cut_mats <- Gam_cut_mats <- vector("list", length = n_components) 

  for (k in seq_len(n_components)) {
    centered_X_F_gam_k <- center_X_F_by_gam(X, Fmat, gams[, k], center_F_marg = center_F_marg)
    X_c_mats[[k]] <- centered_X_F_gam_k$X_c
    Fmat_c_mats[[k]] <- centered_X_F_gam_k$Fmat_c
    
    # # get prepared for computing gam using selected variables
    # sele_var_k <- get_ind_row(Gam_mats[[k]])
    # X_cut_mats[[k]] <- X_c_mats[[k]][, sele_var_k, drop = FALSE] # n * p_sele
    # Gam_cut_mats[[k]] <- Gam_mats[[k]][sele_var_k, , drop = FALSE] # p_sele * d
    
  }

  
  # compute gams 
  for (i in seq_len(n)) {
    num_deno <- rep(0, n_components) # record the number in the denominator
    for (j in seq_len(n_components)) {
      for (k in seq_len(n_components)) {
        if (k == j) {
          temp <- 0
        } else {
          temp <- X_c_mats[[j]][i,,drop=FALSE] - 0.5*matrix(Fmat_c_mats[[j]][i, ], nrow = 1) %*% 
                           t(Gam_mats[[k]]%*%eta_mats[[k]] + Gam_mats[[j]]%*%eta_mats[[j]])
          temp <- matrix(temp, nrow = 1)
          temp <- 1/sigma2 * temp %*% (Gam_mats[[k]]%*%eta_mats[[k]] - Gam_mats[[j]]%*%eta_mats[[j]]) %*% 
                  matrix(Fmat_c_mats[[j]][i, ], ncol = 1)
        }
        num_deno[j] <- num_deno[j] + weight_pi[j] * exp(temp)
      }
      gams[i, j] <- weight_pi[j] / num_deno[j]
    }
  }
  
  # if NA in gams
  if (sum(is.na(gams)) != 0) {
    cat("E step, update gams, number of NA: ", sum(is.na(gams)), "\n")
    gams[is.na(gams)] <- 0
  }
  
  # ----- compute gams, only using the selected variables ----- #
  # gams_sele <- gams
  # for (i in seq_len(n)) {
  #   num_deno <- rep(0, n_components) # record the number in the denominator
  #   for (j in seq_len(n_components)) {
  #     for (k in seq_len(n_components)) {
  #       if (k == j) {
  #         temp <- 0
  #       } else {
  #         temp <- X_cut_mats[[j]][i,,drop=FALSE] - 0.5*matrix(Fmat_c_mats[[j]][i, ], nrow = 1) %*% 
  #           t(Gam_cut_mats[[k]]%*%eta_mats[[k]] + Gam_cut_mats[[j]]%*%eta_mats[[j]])
  #         temp <- matrix(temp, nrow = 1)
  #         temp <- 1/sigma2 * temp %*% (Gam_cut_mats[[k]]%*%eta_mats[[k]] - Gam_cut_mats[[j]]%*%eta_mats[[j]]) %*% 
  #           matrix(Fmat_c_mats[[j]][i, ], ncol = 1)
  #       }
  #       num_deno[j] <- num_deno[j] + weight_pi[j] * exp(temp)
  #     }
  #     gams_sele[i, j] <- weight_pi[j] / num_deno[j]
  #   }
  # }
  # 
  # # if NA in gams
  # if (sum(is.na(gams_sele)) != 0) {
  #   cat("E step, update gams, number of NA: ", sum(is.na(gams)), "\n")
  #   gams_sele[is.na(gams_sele)] <- 0
  # }
  
  # -----               compare the difference            ----- #
  
  
  ### compute log-likelihood
  loglk <- 0
  if (compute_loglk) {
    # under high dimensions, it is computational expensive to compute the 
    # likelihood. Thus, we compute log-likelihood using selected variables
    pi_hat <- colMeans(gams)
    
    for (i in seq_len(n)) {
      tmp <- 0
      for (j in seq_len(n_components)) {
        sele_var <- get_ind_row(Gam_mats[[j]])
        p_sele <- length(sele_var)
        sigma_hat <- sigma2 * diag(p_sele)
        
        temp_mean <- Fmat_c_mats[[j]][i, ,drop = FALSE] %*% t(Gam_mats[[j]]%*%eta_mats[[j]]) # 1 * p
        logNj <- compute_logdmvnorm(X_c_mats[[j]][i, sele_var, drop = FALSE],
                                    temp_mean[sele_var],
                                    sigma_hat)
        tmp <- tmp + pi_hat[j] * exp(logNj)
      }
      loglk <- loglk + log(tmp)
    }
  }
  ###
  
  list(gams = gams, loglk = loglk)
}

# M-step: Update the parameters based on the expected values
m_step_iso_K <- function(X, Y, gams, lam_factor, plot, FUN, 
                         d_vec = d_vec, center_F_marg = FALSE) {
  # Implement the M-step to update parameters based on expected values of latent variables
  # X: n*p data matrix
  # Y: n dimensional vector
  # gams: n*n_components matrix
  # lam_factor: a vector consists of lam_factor of each component
  #             each lam will pass into sparsepca::spca function
  #             as the value of the `alpha` parameter
  # d_vec: a vector of true dimension of central subspaces
  # plot: boolean, if plot dcor plot in the turning process
  # FUN: fitting function
  
  # update estimations of subspaces
  Fmat <- compute_Fmat(Y, FUN) # get Fmat
  n <- nrow(X)
  p <- ncol(X)
  n_componets <- length(d_vec) # get the number of clusters
  D_mats <- vector("list", n_componets)
  for (k in seq_len(n_componets)) {
    D_mats[[k]] <- diag(gams[, k])
  }

  
  # compute matrices before updating
  Sig_fit_mats <- x_c_mats <- Fmat_c_mats <- vector("list", n_componets)
  Gam_mats <- eta_mats <- vector("list", n_componets)
  Sig_mats <- P_Gam_mats <- vector("list", n_componets) # the sample covariance matrix of each cluster
  
  
  for (k in seq_len(n_componets)) {
    MU_out <- MU(X, Y, type = 'mixpfc', FUN = FUN, gam = gams[, k], 
                  center_F_marg = center_F_marg, sig_fit = TRUE)
    Sig_fit_mats[[k]] <- MU_out$sigma_fit
    Sig_mats[[k]] <- MU_out$M
    x_c_mats[[k]] <- MU_out$x_c
    Fmat_c_mats[[k]] <- MU_out$Fmat_c
  }
  
  if (all(lam_factor != 0)) {
    # use sparsepca::spca function
    for (k in seq_len(n_componets)) {
      Gam_mats[[k]] <- sparsepca::spca(Sig_fit_mats[[k]], k = d_vec[k], 
                                       verbose = FALSE, alpha = lam_factor[k], 
                                       tol = 1e-4)$loadings
    }
   
  } else if (all(lam_factor == 0)) {
    # update gam
    for (k in seq_len(n_componets)) {
      Gam_mats[[k]] <- eigen(Sig_fit_mats[[k]], symmetric = TRUE)$vectors[, 1:d_vec[k], drop = FALSE]
    }
  } 
  
  
  
  # update eta
  for (k in seq_len(n_componets)) {
    eta_mats[[k]] <- t(Gam_mats[[k]]) %*% t(x_c_mats[[k]]) %*% D_mats[[k]] %*% 
                     Fmat_c_mats[[k]] %*% ginv(t(Fmat_c_mats[[k]]) %*% D_mats[[k]] %*% Fmat_c_mats[[k]])
  }
  
  
  # update sigma^2
  for (k in seq_len(n_componets)) {
    P_Gam_mats[[k]] <- Gam_mats[[k]] %*% t(Gam_mats[[k]])
  }
  
  sigma2 <- 0
  for (k in seq_len(n_componets)) {
    sigma2 <- sigma2 + tr(D_mats[[k]]) * (tr(Sig_mats[[k]]) - tr(P_Gam_mats[[k]]%*%Sig_fit_mats[[k]]))
  }
  sigma2 <- 1 / (n*p) * sigma2
  
  
  # update weights
  weight_pi <- colMeans(gams)
  
  list(Gam_mats = Gam_mats, weight_pi = weight_pi, gams = gams,
       sigma2 = sigma2, eta_mats = eta_mats)
}




# EM algorithm
mixPFC_iso_K <- function(X, Y, d_vec, FUN = NULL, 
                         n_components = length(d_vec), gams_init = NULL, 
                         lam_factor = rep(0.003, n_components), 
                         weight_init = rep(1/n_components, n_components),
                         center_F_marg = FALSE, compute_loglk = FALSE,
                         max_iter = 100, tol = 1e-3,
                         component_ind = NULL, compute_dc = FALSE,
                         plot = FALSE, print_info = FALSE) {
  # === Input ===
  # X: n*p data matrix
  # Y: n-dimensional observation vector for response
  # d_vec: a vector of true dimensions
  # lam_factor: parameter that pass to sparsepca::spca
  # other parameters have the same meaning as mixPFC_K function
  n <- nrow(X) # number of samples
  p <- ncol(X) # number of predictors
  
  Fmat <- compute_Fmat(Y, FUN) # get non-centered Fmat
 
  
  ## initialization
  weight_pi <- weight_init # Initialize weights 
  Gam_mats <- vector("list", n_components)
  for (i in seq_len(n_components)) {
    Gam_mats[[i]] <- matrix(0, p, d_vec[i]) # Initialize Gam matrices
  }

  # initialize gam
  if (is.null(gams_init)) {
    gams <- matrix(0, n, n_components)
    for (i in seq_len(n)) {
      random_col <- sample(seq_len(n_components), size = 1)
      gams[i, random_col] <- 1
    }
  } else {
    gams <- gams_init
  }
  
  
  # iteration starts
  iteration <- 1
  converged <- FALSE
  loglk <- rep(NA, max_iter)
  distance_correlation <- matrix(NA, nrow = max_iter, ncol = n_components) # distance correlation
  er <- rep(NA, max_iter) # vector to store error rate
  while (iteration <= max_iter && !converged) {
    
    # M-step
    updated_params <- m_step_iso_K(X, Y, gams, lam_factor = lam_factor,
                                   plot = plot, FUN = FUN, d_vec = d_vec,
                                   center_F_marg = center_F_marg)
    # compute the distance correlation between Y and current estimated subspace
    # record these to see if there's improvement 
    if (compute_dc) {
      for (k in seq_len(n_components)) {
        distance_correlation[iteration, k] <- 
          dcor(diag(sqrt(gams[, k])) %*% X %*% updated_params$Gam_mats[[k]], 
               diag(sqrt(gams[, k])) %*% matrix(Y, ncol = 1))
      }
    }
    
    # record error rate 
    if (!is.null(component_ind)) {
      comp_hat <- est_component(gams)
      er[iteration] <- error_rate_K(comp_hat, component_ind, n_components)$error_rate
    }
    
    
    # E-step
    gams_loglk <- e_step_iso_K(X, Fmat, updated_params, 
                               center_F_marg = center_F_marg,
                               compute_loglk = compute_loglk)
    gams <- gams_loglk$gams
    loglk[iteration] <- gams_loglk$loglk
    
    ### compute the difference between current and previous parameters
    diff_Gam_mat <- rep(NA, n_components) 
    diff_weights <- NA
    for (k in seq_len(n_components)) {
      diff_Gam_mat[k] <- max(abs(updated_params$Gam_mats[[k]] - Gam_mats[[k]]))
      diff_weights <- max(abs(updated_params$weight_pi - weight_pi))
    }
    
    ######### print 
    #cat("dims of updated_params$B1: ", dim(updated_params$B1), "\n")
    #cat("dims of B1: ", dim(B1), "\n")
    if (print_info) {
      cat("iteration:", iteration, "\n")
      cat("sigma2:", updated_params$sigma2, "\n")
      for (k in seq_len(n_components)) {
        cat(paste0("difference Gam", k, ":"), diff_Gam_mat[k], " ")
      }
      cat("\n")
      cat("difference weights:", diff_weights, "\n")
      cat("loglk:", loglk[iteration], "\n")
    }
    #########
    ## Check for convergence
    if (all(diff_Gam_mat < tol) && diff_weights < tol) {
      converged <- TRUE
    }

    # if not converge, update params
    Gam_mats <- updated_params$Gam_mats
    weight_pi <- updated_params$weight_pi
    sigma2 <- updated_params$sigma2
    eta_mats <- updated_params$eta_mats
    
    iteration <- iteration + 1
  }
  
  ### compute cluster means 
  mu <- compute_mu_K(X, gams)
  ###
  
  
  list(Gam_mats = Gam_mats, weight_pi = weight_pi, gams = gams, loglk = loglk,
       distance_correlation = distance_correlation, er = er, sigma2 = sigma2, 
       eta_mats = eta_mats, mu = mu)
}


