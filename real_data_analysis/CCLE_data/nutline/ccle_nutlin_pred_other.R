source("LassoSIR_revised.R")
source("seas.R")
source("utility.R")

library(mnormt)
library(energy)
library(msda)
library(mvtnorm)
library(here)
library(glmnet)

# load data
load("/ccle_nutlin.Rdata")

# load test index
# first use gen_index_nutlin.R to generate the index file
load("test_ind_mat_nutlin.Rdata")

N <- nrow(X) # number of samples
p <- ncol(X) # number of variables
ratio <- 0.2
ntest <- floor(N * ratio) # use 1/5 as test data
ntrain <- N - ntest
FUN <- function(y) c(y, y^2) # or y^3
d_subspace <- 1
H <- 5

source("foo2.R")

library(foreach)
library(doParallel)
registerDoParallel(25)
RNGkind("L'Ecuyer-CMRG")
set.seed(2023)
n_sim <- 100
pred_results <- foreach(i=1:n_sim) %dopar% foo2(i) 

file_name <- paste0(paste("nutlin", "other", "pred", sep = "_"), "_y^2.Rdata")

save(pred_results, file = paste0(here(), "/", "results/", file_name))