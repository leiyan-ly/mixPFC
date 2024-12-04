scr_folder <- "mixPFC_src/"
library(mnormt)
library(energy)
library(msda)
library(mvtnorm)
library(here)

source(paste0(scr_folder, "seas.R"))
source(paste0(scr_folder, "utility.R"))
source(paste0(scr_folder, "mixPFC_K.R"))


# load data
load("ccle_nutlin.Rdata")

# load test index
# first use gen_index_nutlin.R to generate the index file
load("test_ind_mat_nutlin.Rdata")


N <- nrow(X) # number of samples
p <- ncol(X) # number of variables
ratio <- 0.2
ntest <- floor(N * ratio) # use 1/5 as test data
ntrain <- N - ntest

K <- 3 # number of clusters
d_dcor <- 2 * floor(ntrain / log(ntrain))
d_pca <- 10
FUN <- function(y) c(y, y^2) 
lam_factor <- rep(0.1, K) 
em_maxiter <- 50
tol_init <- 10e-3
d_subspace <- 1


source("foo1.R")

library(foreach)
library(doParallel)
registerDoParallel(25)
RNGkind("L'Ecuyer-CMRG")
set.seed(2023)
n_sim <- 100
pred_results <- foreach(i=1:n_sim) %dopar% foo1(i) 

file_name <- paste0(paste("nutlin", "K", K, "mixpfc", "pred", sep = "_"), "_y^2.Rdata")

save(pred_results, file = paste0(here(), "/", "results/", file_name))