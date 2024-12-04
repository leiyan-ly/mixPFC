## Main Functions

The implementation of the proposed `mixPFC` algorithm is available in
the folder **mixPFC_src**. The primary functions are contained in the R
files `mixPFC_K.R` and `tune_K.R`, which implement the core `mixPFC`
algorithm and the method for selecting the number of clusters (*K*),
respectively.

### Additional Resources:

-   **Utility Functions**: Additional R files in **mixPFC_src** provide
    utility functions and the implementation of the `mixPFC-ISO`
    algorithm.  
-   **Simulations and Data Analysis**:
    -   Code for simulations of model M4 is provided in the folder
        **simulation**.  
    -   Scripts for analyzing two real datasets are included in
        **real_data_analysis**.

The sections below provide a detailed description of the key functions
`mixPFC_K()` and `tune_K_gap()`.

### `mixPFC_K()`

The `mixPFC_K()` function estimates the parameters of the `mixPFC`
model.

#### Arguments

-   **`X`**: An *n* × *p* matrix of predictor variables.
-   **`Y`**: A length-*n* vector of univariate responses.
-   **`n_components`**: The number of clusters, *K*.
-   **`foldid`**: A vector of integers identifying the fold assignment
    for each observation.
-   **`FUN`**: A *q*-dimensional fitting function, **f**.
-   **`gams_init`**: Initial values for the sample membership weights,
    *γ*<sub>*i**w*</sub>.
-   **`weight_init`**: Initial values for the cluster weights,
    *π*<sub>*w*</sub>.
-   **`lam_factor`**: The tuning parameter vector.
-   **`max_iter`**: Maximum number of iterations allowed.
-   **`tol`**: Stopping threshold for the `mixPFC` algorithm.
-   **`print_info`**: If `TRUE`, prints the error at each iteration.
-   **`compute_dc`**: If `TRUE`, computes the distance correlation at
    each iteration.
-   **`component_ind`**: A vector of true cluster labels (default is
    `NULL`).
-   **`plot`**: If `TRUE`, plots the evaluation results for each tuning
    parameter in the `msda` function.
-   **`compute_loglk`**: If `TRUE`, computes the log-likelihood at each
    iteration.
-   **`perturb`**: If the estimated covariance matrix **Δ**<sup>\*</sup>
    is rank-deficient, adds perturb × **I** to **Δ**<sup>\*</sup> for
    log-likelihood computation. Used only when `compute_loglk = TRUE`.

#### Value

-   **`B_mats`**: A length-*K* list of estimated **B**<sub>*w*</sub>
    matrices.
-   **`weight_pi`**: A length-*K* vector of estimated weights
    *π*<sub>*w*</sub>.
-   **`gams`**: An *n* × *K* matrix of membership weights for each
    sample.
-   **`A_mats`**: A length-*K* list of estimated **A**<sub>*w*</sub>
    matrices.
-   **`Delta_inv`**: An *q**K* × *q**K* matrix of the estimated inverse
    **Δ**<sup>\*</sup>.
-   **`mu`**: A length-*K* list of estimated **μ**<sub>*w*</sub>.
-   **`loglk`**: A length-`max_iter` vector of log-likelihood values.
-   **`distance_correlation`**: A `max_iter`-by-*K* matrix of distance
    correlations.
-   **`er`**: A length-`max_iter` vector of error rates (not `NA` if
    `component_ind` is not `NULL`).

------------------------------------------------------------------------

### `tune_K_gap()`

The `tune_K_gap()` function selects the number of clusters, *K*, using
gap statistics.

#### Arguments

-   **`K_max`**: Maximum number of clusters to consider.
-   **`X`**: An *n* × *p* matrix of predictor variables.
-   **`Y`**: A length-*n* vector of univariate responses.
-   **`d`**: The number of dimensions for each subspace (assumes all
    clusters have the same dimension).
-   **`lam_factor`**: A tuning parameter.
-   **`n_B`**: Number of random samples used to compute the gap
    statistics. If `n_B = 0`, the within-cluster dispersion
    *V*<sub>*K*</sub> is used to select *K*.
-   **`within_clus_metric`**: Metric used for within-cluster evaluation.
    Options include:
    -   `ss_proj_orth`: Sum of squares of projected predictors
        **Q**<sub>*w*</sub>**X**.
    -   `angle`: Sum of angles.
    -   `ss_xy`: Sum of squares of **X** ∘ *Y*.
    -   `ss`: Sum of squares.
-   **`ref_dist`**: Distribution for sampling reference data. Options:
    `scaledPCA`, `original`.
-   **`d_dcor`**: Number of variables to retain during distance
    correlation screening.
-   **`d_pca`**: Number of principal components.
-   **`tol_init`**: Early stopping threshold.
-   **`component_ind`**: Vector of true cluster labels.
-   **`true_comp_ind`**: If `TRUE`, uses `component_ind` to generate
    initial values.
-   **`max_iter`**: Maximum number of iterations allowed.
-   **`tol`**: Stopping threshold for the `mixPFC` algorithm.

#### Value

-   **`gap_stat`**: A list with the following components:
    -   `Tab`: A matrix with `K_max` rows and four columns: `"logW"`,
        `"E.logW"`, `"gap"`, and `"SE.sim"`, where `gap = E.logW - logW`
        and `SE.sim` is the standard error of the gap.
    -   `spaceH0`: The value of the `ref_dist` argument.
    -   `n`: Number of samples.
    -   `B`: Number of random datasets sampled from `ref_dist`.
    -   `call`: The `clusGap_mixpfc()` call.
-   **`er_init`**: A vector of initial error rates (not `NA` if
    `component_ind` is not `NULL`).

## Simulation Illustration of `mixPFC_K()`

This section demonstrates how to reproduce a single replicate of results
under the multiple-cluster model M4.

### Steps to Reproduce:

1.  **Load Required Packages and Source Files**  
    Begin by loading the necessary R packages and sourcing the required
    functions.

``` r
rm(list = ls())
library(mnormt)
library(energy)
library(msda)
library(mvtnorm)
library(R.matlab)
source("seas.R")
source("utility.R")
source("mixPFC_K.R")
source("simulation_data.R")
set.seed(2023)
```

1.  **Set Model Parameters and Generate Data**  
    Define the model parameters and generate data with *K* = 3 clusters
    and a covariance matrix of AR(0.3).  
    To experiment with other covariance structures, such as 0.1**I**,
    **I**, or AR(0.5), simply adjust the value of the variable `j`.

``` r
s <- 10 # number of important variables
K <- 3  # number clusters
N <- 200 * K # sample size
p <- 1000 # number of predictors
FUN <- function(y) c(y, y^2, y^3) # fitting function
beta <- beta_gen_K_random(s, K, 2, seed_ = 2023, print_dist = TRUE)
beta <- lapply(beta, add_zeros, p-s)
ind_row <- lapply(beta, get_ind_row)
ind_row_union <- Reduce(union, ind_row) # the union of important variables of all cluster
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
lam_factor <- list(rep(0.002, K), 
                   rep(0.04, K),
                   rep(0.012, K),
                   rep(0.012, K))
weight <- rep(1 / K, K) # weights vector pi

# generate data with covariance structure AR(0.3, p)
j <- 3 # change this value to generate data with other covariance 
Delta <- Delta_list[[j]]
Gamma <- Gamma_list[[j]]
train <- data_gen_K_random(N, p, Delta, Gamma, weight, 
                           center_F_marg = TRUE, K = K)
x_train <- train$x
y_train <- train$y
component_ind <- train$component_ind
```

1.  **Generate Initial Values**  
    Transform the data into a lower-dimensional space using:
    -   **Distance Correlation (dcor)**, and  
    -   **Principal Component Analysis (PCA)**. Then apply the `mixPFC`
        algorithm with early stopping to obtain initial parameter
        estimates.

``` r
d_dcor <- rep(2*floor(N / log(N)), n_Delta)
d_pca <- rep(10, n_Delta)
tol_init <- c(12e-3, 12e-3, 12e-3, 12e-3)
er_init <- er_kmeanX_reduced_Y <- rep(NA, n_Delta)
em_maxiter <- 40 # maximum number of iteration of mixPFC, this is not the number of short runs that used in initialization
er_rate <- rep(NA, n_Delta)
pi_hat <- matrix(nrow = n_Delta, ncol = K)
TPR <- FPR <- matrix(nrow = n_Delta, ncol = K)
TPR_union <- FPR_union <- rep(NA, n_Delta)
subspace_dist <- matrix(nrow = n_Delta, ncol = K)

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
```

1.  **Run the `mixPFC` Algorithm**  
    With the data prepared and initial values obtained, proceed to run
    the `mixPFC` algorithm.

``` r
nfolds <- 5
foldid <- sample(rep(seq(nfolds), length.out = N))
# run mixpfc with gams_init
mixpfc_fit <- mixPFC_K(x_train, y_train, foldid, FUN = FUN,
                       n_components = K, gams_init = gams_init,
                       lam_factor = lam_factor[[j]], weight_init = colMeans(gams_init),
                       center_F_marg = TRUE, compute_loglk = FALSE,
                       print_info = FALSE, max_iter = em_maxiter)

# eval results
component_ind_hat <- est_component(mixpfc_fit$gams)


er_and_ord <- error_rate_K(component_ind_hat, component_ind, K) # error rate and the best permutation of cluster
# store the error rates
er_rate[j] <- er_and_ord$error_rate
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
```

## Simulation Illustration of `tune_K_gap()`

This section demonstrates how to select the number of clusters using the
method proposed in the paper.

### Steps to Reproduce:

1.  **Load Required Packages and Source Files**  
    Begin by loading the necessary R packages and sourcing the required
    functions.

``` r
library(mnormt)
library(energy)
library(msda)
library(mvtnorm)
source(paste0(scr_folder, "seas.R"))
source(paste0(scr_folder, "utility.R"))
source(paste0(scr_folder, "mixPFC_K.R"))
source(paste0(scr_folder, "simulation_data.R"))
source(paste0(scr_folder, "tune_K.R"))
set.seed(2023)
```

1.  **Generate Data** and **Set Parameters for Cluster Selection**  
    Simulate data under model M4 with *K* = 3 clusters and a covariance
    matrix of AR(0.3).  
    To experiment with other covariance structures—such as 0.1**I**,
    **I**, or AR(0.5)—adjust the value of the variable `j`.  
    Specify the maximum number of clusters to evaluate, setting it to 10
    in this example.

``` r
s <- 10 
K <- 3 # number clusters
N <- 200 * K
p <- 1000
FUN <- function(y) c(y, y^2, y^3) # change this value
# beta <- beta_gen_K(s, sim_set, K, scenario)
beta <- beta_gen_K_random(s, K, 2, seed_ = 2023, print_dist = TRUE)
beta <- lapply(beta, add_zeros, p-s)
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


lam_factor <- c(0.002, 0.04, 0.012, 0.012)
weight <- rep(1 / K, K)
d_dcor <- rep(2*floor(N / log(N)), n_Delta)
d_pca <- rep(10, n_Delta)
tol_init <- c(12e-3, 12e-3, 12e-3, 12e-3)
em_maxiter <- 40

# generate data with covariance structure AR(0.3, p)
j <- 3 # change this value to generate data with other covariance 
Delta <- Delta_list[[j]]
Gamma <- Gamma_list[[j]]
train <- data_gen_K_random(N, p, Delta, Gamma, weight, 
                           center_F_marg = TRUE, K = K)
x_train <- train$x
y_train <- train$y
component_ind <- train$component_ind
```

1.  **Run the `tune_K_gap()` Function**  
    Use the `tune_K_gap()` function to compute gap statistics and
    identify the optimal number of clusters.

``` r
K_max <- 10
tune_result <- tune_K_gap(K_max, d_vec[1], x_train, y_train, component_ind, lam_factor[j], 
                          n_B = 0, within_clus_metric = "ss_proj_orth", true_comp_ind = FALSE,
                          ref_dist = "original", max_iter = em_maxiter,
                          d_dcor = d_dcor[j], d_pca = d_pca[j], tol_init = tol_init[j])
```

By following these steps, you can replicate the process of selecting the
number of clusters and assess the performance of the `tune_K_gap()`
function under different simulation scenarios.
