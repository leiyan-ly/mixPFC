scr_folder <- "mixPFC_src"
library(mnormt)
library(energy)
library(msda)
library(mvtnorm)
library(R.matlab)

source(paste0(scr_folder, "utility.R"))
source(paste0(scr_folder, "simulation_data.R"))
set.seed(2023)


#### code for generate data with 3 clusters ####
save_data_path <- "save_data/3_cluster/"
cov_path <- c("cov_01I", "cov_I", "cov_AR03", "cov_AR05")

s <- 10
K <- 3 # number clusters
N <- 200 * K
p <- 1000
FUN <- function(y) c(y, y^2, y^3) # change this value
# beta <- beta_gen_K(s, sim_set, K, scenario)
beta <- beta_gen_K_random(s, K, 2, seed_ = 2023, print_dist = TRUE)
beta <- lapply(beta, add_zeros, p-s)

Delta_list <- list(0.1*diag(p), diag(p), AR(0.3, p), AR(0.5, p))
n_Delta <- length(Delta_list)
Gamma_list <- vector("list", n_Delta)
for (i in seq_len(n_Delta)) { # create Gamma for each Delta
  Gamma <- vector("list", length = K)
  for (k in seq_len(K)) {
    if (i <= 2) {
      Gamma[[k]] <- beta[[k]]
    } else {
      Gamma[[k]] <- Delta_list[[i]] %*% beta[[k]]
    }
    
  }
  Gamma_list[[i]] <- Gamma
}

weight <- rep(1 / K, K)

############################################################
# use the following foo function to generate and save data #
############################################################
foo <- function(i) {
  set.seed(i)
  for (j in seq_len(n_Delta)) {
    
    Delta <- Delta_list[[j]]
    Gamma <- Gamma_list[[j]]
    train <- data_gen_K_random(N, p, Delta, Gamma, weight, 
                               center_F_marg = TRUE, K = K)
    #### use R.matlab to save the generated data #####
    file_path <- paste0(save_data_path, cov_path[j], '/', 'train_', i, '.mat')
    writeMat(file_path, train = train)
  }
}


library(foreach)
library(doParallel)
registerDoParallel(20)
RNGkind("L'Ecuyer-CMRG")
set.seed(2023)
n_sim <- 100
sim_result <- foreach(i=1:n_sim) %dopar% foo(i)

## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##

#### code for generate data with 5 clusters ####

save_data_path <- "save_data/5_cluster/"
cov_path <- c("cov_01I", "cov_I", "cov_AR03", "cov_AR05")

s <- 10
K <- 5 # number clusters
N <- 200 * K
p <- 1000
FUN <- function(y) c(y, y^2, y^3) # change this value
# beta <- beta_gen_K(s, sim_set, K, scenario)
beta <- beta_gen_K_random(s, K, 2, seed_ = 2023, print_dist = TRUE)
beta <- lapply(beta, add_zeros, p-s)

Delta_list <- list(0.1*diag(p), diag(p), AR(0.3, p), AR(0.5, p))
n_Delta <- length(Delta_list)
Gamma_list <- vector("list", n_Delta)
for (i in seq_len(n_Delta)) { # create Gamma for each Delta
  Gamma <- vector("list", length = K)
  for (k in seq_len(K)) {
    if (i <= 2) {
      Gamma[[k]] <- beta[[k]]
    } else {
      Gamma[[k]] <- Delta_list[[i]] %*% beta[[k]]
    }
    
  }
  Gamma_list[[i]] <- Gamma
}

weight <- rep(1 / K, K)

############################################################
# use the following foo function to generate and save data #
############################################################
foo <- function(i) {
  set.seed(i)
  for (j in seq_len(n_Delta)) {
    
    Delta <- Delta_list[[j]]
    Gamma <- Gamma_list[[j]]
    train <- data_gen_K_random(N, p, Delta, Gamma, weight, 
                               center_F_marg = TRUE, K = K)
    #### use R.matlab to save the generated data #####
    file_path <- paste0(save_data_path, cov_path[j], '/', 'train_', i, '.mat')
    writeMat(file_path, train = train)
  }
}


library(foreach)
library(doParallel)
registerDoParallel(20)
RNGkind("L'Ecuyer-CMRG")
set.seed(2023)
n_sim <- 100
sim_result <- foreach(i=1:n_sim) %dopar% foo(i)

