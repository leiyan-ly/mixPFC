# This file generates data for mixture PFC model
#rm(list = ls())
#library(pbmcapply)
library(MASS)
#source("utility.R")
# #############################################

beta_gen_K_random <- function(p, K, d, seed_ = NULL, print_dist = FALSE) {
  # p: dimension of predictor, an even number 
  # K: number of clusters
  # d: dimension of each subspace
  # print_dist: bool, if the print the subspace distance between subspaces
  if (!is.null(seed_))
    set.seed(seed_)
  # define the range of uniform distribution
  lb <- 2.1
  up <- 2.5
  betaK <- vector("list", K)
  for (k in seq_len(K)) {
    beta_sign <- matrix(sample(c(-1, 1), p*d, replace = TRUE, prob = c(0.5, 0.5)),
                        nrow = p)
    betaK[[k]] <- matrix(runif(p*d, min = lb, max = up), nrow = p) * beta_sign
  }
  
  if (print_dist) {
    subspace_dist <- matrix(nrow = K, ncol = K)
    for (i in seq_len(K-1)) {
      for (j in seq(i+1, K, 1)) {
        subspace_dist[i, j] <- subspace(betaK[[i]], betaK[[j]])
      }
    }
    print(subspace_dist)
  }
  
  betaK
}



data_gen_K_random <- function(N, p, Delta, Gamma, weight, 
                              mu = NULL, center_F_marg = FALSE,
                              K = NULL){                
  # Data generating function, assume all subspaces have same dimension
  # N: number of samples
  # p: number of predictors
  # Delta: p*p covariance matrix
  # Gamma: a list of p*d matrix
  # weight: a vector of cluster weights
  # sim_set: simulation setting, an integer
  # mu: a list of vectors, contains intercept for each PFC
  d_vec <- sapply(Gamma, ncol)
  error <- mvrnorm(N, rep(0, p), Delta)
  x <- matrix(0, N, p)
  component_ind <- rep(0, N)
  # if mu is not provided, set to 0
  if (is.null(mu)) { 
    mu <- vector("list", K)
    for (k in seq_len(K)) {
      mu[[k]] <- rep(0, p)
    }
  }
  
  if (all(d_vec == 1)) {
    # each cluster is one dimension
      y <- runif(N, -1, 1)
      f <- cbind(y, abs(y))
      etas <- vector("list", K)
      for (k in seq_len(K)) {
        etas[[k]] <- matrix(c(1, 0.3), nrow = 1)
      }

      ### center f(y)
      if (center_F_marg) {
        f_mean <- colMeans(f)
        #cat("F mean:", f_mean, "\n")
        f <- f - matrix(f_mean, nrow(f), ncol(f), byrow = TRUE) # centered function f
      }
      ###
      for (i in seq_len(N)) {
        component <- sample(1:K, size = 1, prob = weight)
        x[i, ] <- mu[[component]] + f[i, , drop = FALSE] %*% t(etas[[component]]) %*% t(Gamma[[component]]) + error[i, ]
        component_ind[i] <- component
      }
    } else if (all(d_vec == 2)) { # each cluster is two dimensions
      y <- runif(N, -1, 1)
      f <- cbind(y, abs(y))
      ### center f(y)
      if (center_F_marg) {
        f_mean <- colMeans(f)
        f <- f - matrix(f_mean, nrow(f), ncol(f), byrow = TRUE) # centered function f
      }
      ###
      for (i in seq_len(N)) {
        component <- sample(1:K, size = 1, prob = weight)
        x[i, ] <- mu[[component]] + f[i, , drop = FALSE] %*% t(Gamma[[component]]) + error[i, ]
        component_ind[i] <- component
      }
      
    }
    
  
  list(x = x, y = y, component_ind = component_ind, error = error)
}

beta_gen3 <- function(p, sim_set) {
  # p: dimension of predictor, an even number 
  # sim_set: simulation setting, an integer
  
  if (sim_set == 9) { # d1 = d2 = 1, same space
    beta1 <- matrix(rep(2.3, p), ncol = 1)
    #beta1 <- matrix(c(runif(6, 2, 2.5), rep(0, p-6)), ncol = 1)
    beta2 <- -beta1
  } else if (sim_set == 10) { # d1 = d2 = 1, orthogonal
    beta1 <- matrix(rep(2.3, p), ncol = 1)
    #beta2 <- matrix(qr.Q(qr(beta1), complete = TRUE)[, 2] * sqrt(sum(beta1^2)), ncol = 1)
    # the problem of QR decomposition is the obtained beta2 has loadings that are not close to each other
    beta2 <- matrix(rep(c(1, -1), times = p/2), ncol = 1)
    # sum(beta1*beta2) is close to 0. we modify the last element of beta2 a little bit
    # such that beta1 and beta2 are orthogonal
    beta2[p] <- -sum(beta1[1:(p-1)] * beta2[1:(p-1)]) / beta1[p]
    beta2 <- beta2 / sqrt(sum(beta2^2)) * sqrt(sum(beta1^2))
    subspace(beta1, beta2)
  } else if (sim_set == 13) { # d1 = d2 = 2, angle
    beta1 <- beta2 <- matrix(0, nrow = p, ncol = 2)
    beta1[, 1] <- rep(2.3, p)
    beta1[, 2] <- rep(2.3, p)
    beta1[, 2][seq(2, p, by = 2)] <- rep(-2.3, p/2)
    
    beta2[, 1] <- -beta1[, 1]
    beta2[1, 1] <- beta1[1, 1]
    beta2[, 1] <- beta2[, 1] / sqrt(sum(beta2[,1]^2)) * sqrt(sum(beta1[,1]^2))
    beta2[, 2] <- -beta1[, 2]
    beta2[1, 2] <- beta1[1, 2]
    beta2[, 2] <- beta2[, 2] / sqrt(sum(beta2[,2]^2)) * sqrt(sum(beta1[,2]^2))
    cat("Simulation setting 13, distance between beta1 and beta2:", subspace(beta1, beta2), "\n")
    
  } else if (sim_set == 11) { # d1 = d2 = 1, angle
    if (p == 6) {
      beta1 <- matrix(c(2.1, 2.3, 2.3, 2.1, 2.4, 2.4), ncol = 1)
      beta2 <- -beta1
      beta2[1] <- beta1[1]
      
    } else {
      stop("beta_gen3, p not equal to 6!")
    }
    
    #subspace(beta1, beta2)
  }
  #cat("beta1: ", beta1, "\n")
  #cat("beta2: ", beta2, "\n")
  cat("D(beta1, beta2): ", subspace(beta1, beta2), "\n")
  
  list(beta1 = beta1, beta2 = beta2)
}

data_gen2 <- function(N, p, Delta, Gamma, weight, mu = list(rep(0, p), rep(0, p)),
                      center_F_marg = FALSE,
                      sim_set = 9, n_components = 2){                # Data generating function
  # N: number of samples
  # sim_set: simulation setting, an integer
  # mu: a 2D vector, contains intercept for each PFC
  
  error <- mvrnorm(N, rep(0, p), Delta)
  x <- matrix(0, N, p)
  component_ind <- rep(0, N)
  if (sim_set == 9 | sim_set == 10 | sim_set == 11) {
    y <- runif(N, -1, 1)
    f <- cbind(y, abs(y))
    eta1 <- eta2 <- matrix(c(1, 0.3), nrow = 1)
    # y <- runif(N, 0, 2)
    # f <- cbind(y^2, y^3)
    # eta1 <- matrix(c(1, 0), nrow = 1)
    # eta2 <- matrix(c(0, 1), nrow = 1)
    
    
    etas <- list(eta1, eta2)
    ### center f(y)
    if (center_F_marg) {
      f_mean <- colMeans(f)
      #cat("F mean:", f_mean, "\n")
      f <- f - matrix(f_mean, nrow(f), ncol(f), byrow = TRUE) # centered function f
    }
    ###
    for (i in seq_len(N)) {
      component <- sample(1:n_components, size = 1, prob = weight)
      x[i, ] <- mu[[component]] + f[i, , drop = FALSE] %*% t(etas[[component]]) %*% t(Gamma[[component]]) + error[i, ]
      component_ind[i] <- component
    }
    
  } else if (sim_set == 12 | sim_set == 13 | sim_set == 1) {
    y <- runif(N, -1, 1)
    f <- cbind(y, abs(y))
    ### center f(y)
    if (center_F_marg) {
      f_mean <- colMeans(f)
      f <- f - matrix(f_mean, nrow(f), ncol(f), byrow = TRUE) # centered function f
    }
    ###
    for (i in seq_len(N)) {
      component <- sample(1:n_components, size = 1, prob = weight)
      x[i, ] <- mu[[component]] + f[i, , drop = FALSE] %*% t(Gamma[[component]]) + error[i, ]
      component_ind[i] <- component
    }
  } else if (sim_set == 14 | sim_set == 15) {
    y <- runif(N, -1, 1)
    f <- cbind(y, abs(y))
    eta2 <- matrix(c(1, 0.3), nrow = 1)
    ### center f(y)
    if (center_F_marg) {
      f_mean <- colMeans(f)
      f <- f - matrix(f_mean, nrow(f), ncol(f), byrow = TRUE) # centered function f
    }
    ###
    for (i in seq_len(N)) {
      component <- sample(1:n_components, size = 1, prob = weight)
      if (component == 1) {
        x[i, ] <- mu[[component]] + f[i, , drop = FALSE] %*% t(Gamma[[component]]) + error[i, ]
      } else {
        x[i, ] <- mu[[component]] + f[i, , drop = FALSE] %*% t(eta2) %*% t(Gamma[[component]]) + error[i, ]
      }
      component_ind[i] <- component
    }
    return(list(x = x, y = y, component_ind = component_ind, error = error))
  }
  
  
  list(x = x, y = y, component_ind = component_ind, error = error)
}
