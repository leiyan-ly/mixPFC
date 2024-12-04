# compute the posterior probability
compute_post_K <- function(X, Y, Beta_mats, weight_pi, gams, FUN = NULL,
                           d_vec = NULL, act_ind = NULL, center_F_marg = FALSE) {
  # after the em algorithm finished, we truncate matrices B1 and B2
  # and then do one more e-step to obtain the posterior probability
  # Inputs:
  # ======
  # X: n*p design matrix
  # Y: length n vector
  # Beta_mats: p*q matrix
  # gams: n*n_components
  # d_vec: true dimensions of Gamma1, Gamma2, ..., GammaK
  # Outputs:
  # =======
  # gams: updated responsibilities
  if (is.null(d_vec)) {
    stop("true dimension d_vec is missing")
  }
  
  Fmat <- compute_Fmat(Y, FUN) 
  q <- ncol(Fmat)
  if (q == nrow(X))
    Fmat <- t(Fmat)
  q <- ncol(Fmat)
  n_components <- length(d_vec)
  beta_hat_mats <- vector("list", length = n_components)
  for (k in seq_len(n_components)) {
    d_k <- d_vec[k]
    beta_hat_mats[[k]] <- svd(Beta_mats[[k]])$u[, 1:d_k, drop = FALSE]
  }
  
  if (is.null(act_ind)) {
    params <- list(B_mats = beta_hat_mats, gams = gams, weight_pi = weight_pi)
    post_prob <- e_step_K(X, Fmat, params, center_F_marg = center_F_marg)$gams
  } else {
    beta_hat_mats_act_ind <- vector("list", length = n_components)
    for (k in seq_len(n_components)) {
      beta_hat_mats_act_ind[[k]] <- beta_hat_mats[[k]][act_ind, , drop = FALSE]
    }
    params <- list(B_mats = beta_hat_mats_act_ind, 
                   gams = gams, weight_pi = weight_pi)
    post_prob <- e_step_K(X[, act_ind], Fmat, params, center_F_marg = center_F_marg)$gams
  }
  
  post_prob
}

# E-step: Compute the expected value of the latent variables
e_step_K <- function(X, Fmat, params, perturb = 1e-4, center_F_marg = FALSE,
                     compute_loglk = FALSE) {
  # Implement the E-step to compute expected values of latent variables
  # Inputs:
  # ======
  # X: n*p design matrix
  # Fmat: F matrix of dimension n*q
  # params: a list of params, consists of 
  #        B_mats: a list of p*q matrix, length eqauls to number of clusters
  #        gams: n*n_components
  #        weight_pi: length n_components vector
  # Outputs:
  # =======
  # gams: updated responsibilities
  # loglk: log likelihood of current estimates
  
  B_mats <- params$B_mats
  gams <- params$gams
  weight_pi <- params$weight_pi
  n_components <- ncol(gams) # the number of clusters
  n <- nrow(gams) # number of samples
  
  Z_mats <- W_mats <- Fmat_c <- vector("list", length = n_components)
  
  for (k in seq_len(n_components)) {
    W_mats[[k]] <- diag(gams[, k])
    centered_X_F_gam_k <- center_X_F_by_gam(X, Fmat, gams[, k], center_F_marg = center_F_marg)
    #cat(typeof(B_mats[[k]]), "\n")
    Z_mats[[k]] <- centered_X_F_gam_k$X_c %*% B_mats[[k]] # n*q
    Fmat_c[[k]] <- centered_X_F_gam_k$Fmat_c
  }
  
  B <- do.call(cbind, B_mats) # p*qK
  Z <- do.call(cbind, Z_mats) # n*qK
  
  A_mats <- Delta_mats <- vector("list", length = n_components)
  for (k in seq_len(n_components)) {
    A_mats[[k]] <- t(Z) %*% W_mats[[k]] %*% Fmat_c[[k]] %*% ginv(t(Fmat_c[[k]]) %*% W_mats[[k]] %*% Fmat_c[[k]]) # qK*q
    Delta_mats[[k]] <- (t(Z) - A_mats[[k]]%*%t(Fmat_c[[k]])) %*% W_mats[[k]] %*% (Z - Fmat_c[[k]]%*%t(A_mats[[k]])) * (1/sum(gams[,k]))
  }
  
  Delta <- sum(gams[, 1]) / n * Delta_mats[[1]]
  for (k in 2:n_components) {
    Delta <- Delta + sum(gams[, k]) / n * Delta_mats[[k]]
  }
  
  
  Delta_inv <- ginv(Delta)
  
  # compute gams 
  for (i in seq_len(n)) {
    num_deno <- rep(0, n_components) # record the number in the denominator
    for (j in seq_len(n_components)) {
      for (k in seq_len(n_components)) {
        if (k == j) {
          temp <- 0
        } else {
          temp <- matrix(Z[i,,drop=FALSE] - 0.5*matrix(Fmat_c[[j]][i, ], nrow = 1)%*%t(A_mats[[k]] + A_mats[[j]]), nrow = 1) # 1*qK
          temp <- temp %*% Delta_inv %*% (A_mats[[k]] - A_mats[[j]]) %*% matrix(Fmat_c[[j]][i, ], ncol = 1)
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
  
  # compute log likelihood
  loglk <- 0
  if (compute_loglk) {
    # is Delta invertible?
    min_eigen_val <- min(eigen(Delta, only.values = TRUE)$values)
    if (min_eigen_val > perturb) {
      perturb <- NULL
    }
    pi_hat <- colMeans(gams)
    for (i in seq_len(n)) {
      tmp <- 0
      for (j in seq_len(n_components)) {
        logNj <- compute_logdmvnorm(Z[i, , drop = FALSE],
                                    Fmat_c[[j]][i, ,drop = FALSE] %*% t(A_mats[[j]]),
                                    Delta, perturb = perturb)
        tmp <- tmp + pi_hat[j] * exp(logNj)
      }
      loglk <- loglk + log(tmp)
    }
  }
  
  
  list(gams = gams, loglk = loglk, A_mats = A_mats, Delta_inv = Delta_inv)
}

# M-step: Update the parameters based on the expected values
m_step_K <- function(X, Y, foldid, gams, lam_factor, plot, FUN, center_F_marg = FALSE,
                     print_info = FALSE) {
  # Implement the M-step to update parameters based on expected values of latent variables
  # X: n*p data matrix
  # Y: n dimensional vector
  # gams: n*n_components matrix
  # lam_factor: a vector consists of lam_factor of each component
  # plot: boolean, if plot dcor plot in the turning process
  # FUN: fitting function
  
  # obtain the number of clusters
  n_components <- ncol(gams)
  mixpfc_fit <- vector("list", n_components)
  
  # update estimations of subspaces
  if (all(lam_factor != 0)) {
    for (k in seq_len(n_components)) {
      mixpfc_fit[[k]] <- cv.msda(X, Y, type = 'mixpfc', foldid = foldid, 
                                 FUN = FUN, center_F_marg = center_F_marg,
                                 lambda.factor = lam_factor[k], 
                                 gam = gams[, k], plot = plot)
    }
    
    if (print_info) {
      for (k in seq_len(n_components)) {
        cat("fit", k, "lambda:", mixpfc_fit[[k]]$lam_max, "\n")
      }
    }
    
  } else if (all(lam_factor == 0)) {
    for (k in seq_len(n_components)) {
      if (print_info) {
        cat("m step with lambda = 0 for cluster ", k, "\n")
      }
      MU_out_k <- MU(X, Y, type = 'mixpfc', FUN = FUN, gam = gams[, k], 
                     center_F_marg = center_F_marg)
      M_k <- MU_out_k$M
      U_k <- MU_out_k$U
      mixpfc_fit_k <- list()
      mixpfc_fit_k$beta <- ginv(M_k) %*% U_k
      
      mixpfc_fit[[k]] <- mixpfc_fit_k
    }
    
  } 
  
  # update weights
  weight_pi <- colMeans(gams)
  
  # obtain beta from mixpfc_fit
  B_mats <- vector("list", length = n_components)
  for (k in seq_len(n_components)) {
    B_mats[[k]] <- mixpfc_fit[[k]]$beta
  }
  
  list(B_mats = B_mats, weight_pi = weight_pi, gams = gams)
}



# # compute Fmat
# compute_Fmat <- function(Y, FUN = NULL) {
#   # Y: n-dimensional observation vector for response
#   # FUN: the user-specified function f in SEAS-PFC. The default is f(y) = (y, y^2, y^3)
#   if (is.null(FUN)) Fmat <- cbind(Y, Y^2, Y^3) # the default function
#   else Fmat <- t(sapply(Y, FUN))
#   
#   # n <- nrow(Fmat)
#   # q <- ncol(Fmat)
#   # Fmat_mean <- colMeans(Fmat)
#   # Fmat_c <- Fmat - matrix(Fmat_mean, n, q, byrow = TRUE) # centered function matrix
#   # lb <- apply(Fmat_c, 2, quantile, 0.1)
#   # ub <- apply(Fmat_c, 2, quantile, 0.9)
#   # for (i in seq_len(q)) {
#   #   Fmat_c[, i] <- sapply(Fmat_c[, i], cut_func, lb[i], ub[i]) # cut extreme value
#   # }
#   #if (nrow(Fmat == 1)) Fmat <- t(Fmat)
#   Fmat
# }

# EM algorithm for multiple clusters
mixPFC_K <- function(X, Y, foldid = NULL, FUN = NULL, 
                     n_components = 2, gams_init = NULL, 
                     perturb = 1e-4, lam_factor = rep(0.03, n_components),
                     weight_init = rep(1/n_components, n_components), 
                     center_F_marg = FALSE, compute_loglk = FALSE,
                     compute_dc = FALSE,
                     plot = FALSE, print_info = FALSE,
                     component_ind = NULL, max_iter = 100, tol = 1e-3) {
  # X: n*p data matrix
  # Y: n-dimensional observation vector for response
  # foldid: foldid, use this to do CV
  # FUN: the function to compute Fmat
  # n_components: integer, the number of clusters
  # gams_init: n*n_component matrix, the initial posterior probabilities
  # perturb: the estimated covariance matrix Delta could be rank-deficient
  #          if that happens, we add perturb*I to Delta to compute loglikelihood.
  #          this parameter is used only when compute_loglk = TRUE
  # compute_dc: bool, if we compute distance covaraince 
  # plot: bool, this parameter is passed to cv.msda
  # print_info: bool, passed to m_step_K, should it print the best lambda for each 
  #             cluster
  # component_ind: the true component labels for the samples
  
  n <- nrow(X) # number of samples
  p <- ncol(X) # number of predictors
 
  #X_c <- center_mat(X) # centered X
  Fmat <- compute_Fmat(Y, FUN) # get centered Fmat
  q <- ncol(Fmat)
  if (q == n)
    Fmat <- t(Fmat)
  q <- ncol(Fmat)
  
  ## initialization
  weight_pi <- weight_init # Initialize weights 
  B_mats <- replicate(n_components, matrix(0, nrow = p, ncol = q), simplify = FALSE) # Initialize B matrices
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
    updated_params <- m_step_K(X, Y, foldid, gams, lam_factor = lam_factor,
                               plot = plot, FUN = FUN, center_F_marg = center_F_marg, 
                               print_info = print_info)
    # compute the distance correlation between Y and current estimated subspace
    # record these to see if there's improvement 
    if (compute_dc) {
      for (k in seq_len(n_components)) {
        distance_correlation[iteration, k] <- 
          dcor(diag(sqrt(gams[, k])) %*% X %*% updated_params$B_mats[[k]], diag(sqrt(gams[, k])) %*% Y)
      }
    }
    
    # record error rate 
    if (!is.null(component_ind)) {
      comp_hat <- est_component(gams)
      er[iteration] <- error_rate_K(comp_hat, component_ind, n_components)$error_rate
    }
    
    # E-step
    gams_loglk <- e_step_K(X, Fmat, updated_params, 
                           perturb = perturb, center_F_marg = center_F_marg,
                           compute_loglk = compute_loglk)
    # update gams
    gams <- gams_loglk$gams
    loglk[iteration] <- gams_loglk$loglk

    
    ### compute the difference between current and previous parameters
    diff_Bmat <- rep(NA, n_components) 
    diff_weights <- NA
    for (k in seq_len(n_components)) {
      diff_Bmat[k] <- max(abs(updated_params$B_mats[[k]] - B_mats[[k]]))
      diff_weights <- max(abs(updated_params$weight_pi - weight_pi))
    }
    
    ######### print 
    #cat("dims of updated_params$B1: ", dim(updated_params$B1), "\n")
    #cat("dims of B1: ", dim(B1), "\n")
    if (print_info) {
      cat("iteration:", iteration, "\n")
      for (k in seq_len(n_components)) {
        cat(paste0("difference B", k, ":"), diff_Bmat[k], " ")
      }
      cat("\n")
      cat("difference weights:", diff_weights, "\n")
      cat("loglk:", loglk[iteration], "\n")

    }
    #########
    ## Check for convergence
    if (all(diff_Bmat < tol) && diff_weights < tol) {
      converged <- TRUE
    }

    # if not converge, update params
    B_mats <- updated_params$B_mats
    weight_pi <- updated_params$weight_pi
    
    iteration <- iteration + 1
  }
  
  ### get A and Delta_inv for prediction of new data points 
  A_mats <- gams_loglk$A_mats
  Delta_inv <- gams_loglk$Delta_inv
  ###
  
  ### compute mu for each cluster
  mu <- compute_mu_K(X, gams)
  ###
  
  ### compute log-likelihood using e_step_K after convergence 
  loglk_converged <- e_step_K(X, Fmat, updated_params, 
                              perturb = perturb, center_F_marg = center_F_marg,
                              compute_loglk = TRUE)$loglk
  loglk <- c(loglk, loglk_converged)
  ###
  
  
  list(B_mats = B_mats, weight_pi = weight_pi, gams = gams, loglk = loglk,
       distance_correlation = distance_correlation, er = er, 
       A_mats = A_mats, Delta_inv = Delta_inv, mu = mu)
}

# Example usage





