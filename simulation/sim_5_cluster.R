scr_folder <- "mixPFC_src/"
library(mnormt)
library(energy)
library(msda)
library(mvtnorm)
library(R.matlab)
library(here)
source(paste0(scr_folder, "seas.R"))
source(paste0(scr_folder, "utility.R"))
source(paste0(scr_folder, "mixPFC_K.R"))
source(paste0(scr_folder, "simulation_data.R"))
source(paste0(scr_folder, "init_utility_K.R"))
set.seed(2023)

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
ind_row <- lapply(beta, get_ind_row)
ind_row_union <- Reduce(union, ind_row) # the union of important varaibles of all cluster
d_vec <- sapply(beta, ncol)
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
lam_factor <- list(rep(0.002, K), # 0.005
                   rep(0.04, K),
                   rep(0.012, K),
                   rep(0.012, K))
weight <- rep(1 / K, K)


# variables that used in foo3.R, when generate initial values 
d_dcor <- rep(2*floor(N / log(N)), n_Delta)
d_pca <- rep(10, n_Delta)
tol_init <- c(12e-3, 5e-3, 3e-3, 2e-3)
er_init <- er_kmeanX_reduced_Y <- rep(NA, n_Delta)
em_maxiter <- 40 # maximum number of iteration of mixPFC, this is not the number of short runs that used in initialization
# declare some variables to store the results
er_rate <- er_rate2 <- er_rate_ora <- rep(NA, n_Delta)
pi_hat <- pi_hat2 <- matrix(nrow = n_Delta, ncol = K)
TPR <- FPR <- matrix(nrow = n_Delta, ncol = K)
TPR_union <- FPR_union <- rep(0, n_Delta)
subspace_dist <- subspace_dist_ora <- matrix(nrow = n_Delta, ncol = K)
ind_hat <- vector("list", n_Delta)

source("foo.R")

j_seq <- c(1, 2, 3, 4)

library(foreach)
library(doParallel)
registerDoParallel(25)
RNGkind("L'Ecuyer-CMRG")
set.seed(2023)
n_sim <- 100
sim_result_3_cluster <- foreach(i=1:n_sim) %dopar% foo(i, em_maxiter = em_maxiter) 
file_name <- paste0(paste(paste0(K, "_cluster"), paste0("p", p), paste0("d", d_vec[1]), sep = "_"), "_AR.Rdata")
save(sim_result_3_cluster, file = paste0(here(), "/", "results/", file_name))
# save the estimated labels that can be used by other methods
comp_hat_file_name <- paste0(paste(paste0(K, "_cluster"), paste0("p", p), paste0("d", d_vec[1]), sep = "_"), "_comp_hat", "_AR.Rdata")
save_comp_hat(sim_result_3_cluster, paste0(here(), "/", "results/", comp_hat_file_name))