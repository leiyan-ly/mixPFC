
# generate the test index of AZD6244
RNGkind("L'Ecuyer-CMRG")
set.seed(2023)
n_rep <- 120
N <- 479
ratio <- 0.2
ntest <- floor(N * ratio) # use 1/5 as test data
gen_bs_index <- function(n_rep, ntest) {
  ind_mat <- matrix(nrow = n_rep, ncol = ntest)
  for (i in seq_len(n_rep)) {
    ind_mat[i, ] <- sort(sample(seq_len(N), size = ntest, replace = FALSE))
  }
  ind_mat
}

test_ind_mat <- gen_bs_index(n_rep, ntest)
save(test_ind_mat, file = "test_ind_mat_AZD6244.Rdata")

