
# this script contains function to select the number of clusters K
tune_K_gap <- function(K_max, d, X, Y, component_ind = NULL, lam_fact = NULL, 
                       n_B = 10, within_clus_metric = "angle", ref_dist = "original",
                       true_comp_ind = TRUE, max_iter = 40,
                       tol = 1e-3, d_dcor = NULL, d_pca = NULL, tol_init = NULL) {
  # input:
  # K_max: the maximum number of clusters
  # d: dimension of subspace, assume all clusters have the same dimension
  # n_B: number of sampling datasets
  # ref_dist: original or scaledPCA
  # lam_fact: a number, not a vector (like in `mixPFC_K`)
  # other parameters are defined in function `mixPFC_K`
  
  n <- nrow(X)
  p <- ncol(X)
  K_seq <- c(2:K_max)
  n_K <- length(K_seq) # number of candidate K
  er_init <- rep(NA, n_K)
  
  nfolds <- 5
  foldid <- sample(rep(seq(nfolds), length.out = n))
  
  gams_init_l <- vector("list", n_K)
  lam_fact_l <- vector("list", n_K)
  # generate gams_init for each given K
  for (i in seq_len(n_K)) {
    K <- K_seq[i]
    # cat("K: ", K, "\n")
    # generate initial value
    # if the true labels are provided
    if (true_comp_ind) {
      K_true <- length(unique(component_ind)) # true K
      if (K == K_true) {
        true_percent <- 0.7
      } else {
        true_percent <- 1
      }
      gams_init <- gen_gams_init_K(component_ind, true_percent = true_percent) # n*K^*
      if (K <= K_true) {
        # gams_init should only have K columns
        gams_init <- gams_init[, 1:K, drop = FALSE]
        for (j in seq_len(n)) {
          if (sum(gams_init[j, ]) == 0) {
            random_cluster <- sample(c(1:K), 1)
            gams_init[j, random_cluster] <- 1
          }
          
        }
      } else { # if K > K_true
        gams_init <- cbind(gams_init, matrix(0, nrow = n, ncol = K - K_true))
        n_each_class <- floor(n / K)
        random_cluster <- sample(c(1:n), n_each_class * (K - K_true)) # random label
        random_cluster <- matrix(random_cluster, nrow = n_each_class)
        for (j in seq_len(K - K_true)) {
          gams_init[random_cluster[, j], ] <- 0
          gams_init[random_cluster[, j], j+K_true] <- 1
        }
      }
      
    } else {
      # if the true labels are not provided
      # use dcor + PCA + mixPFC
      if (is.null(d_dcor) | is.null(d_pca) | is.null(tol_init)) {
        stop("missing values: d_dcor, d_pca or tol_init!", "\n")
      }
      dcor_vec <- apply(X, 2, function(x) dcor(x, Y))
      select_ind <- tail(order(dcor_vec), d_dcor)
      x_train_cut <- X[, select_ind, drop = FALSE]
      pca_fit <- prcomp(x_train_cut, scale. = FALSE, rank. = d_pca)
      x_train_cut_reduced <- pca_fit$x
      # kmeans on x_train_cut_reduced*Y
      xy_train <- x_train_cut_reduced * matrix(Y, nrow = nrow(x_train_cut_reduced), ncol = ncol(x_train_cut_reduced))
      kmeanXY_comp_ind <- kmeans(xy_train, K, nstart = 20)$cluster
      gams_init <- gen_gams_init_K(kmeanXY_comp_ind, true_comp_ind = FALSE)
      # short run of mixPFC on the reduced data
      em_init <- mixPFC_K(x_train_cut_reduced, Y, FUN = FUN, 
                          n_components = K, gams_init = gams_init,
                          lam_factor = rep(0, K), weight_init = colMeans(gams_init), 
                          center_F_marg = TRUE, compute_loglk = FALSE,
                          max_iter = 15, tol = tol_init)
      gams_init <- em_init$gams
      
    }
    
    gams_init_l[[i]] <- gams_init
    if (!is.null(component_ind)) {
      er_init[i] <- error_rate_K(est_component(gams_init), component_ind, K)$error_rate
    }
    lam_fact_l[[i]] <- rep(lam_fact, K)
  }
  
  # use gap statistics
  gap_stat <- clusGap_mixpfc(X, Y, K_max, gams_init_l, lam_fact_l, 
                             within_clus_metric = within_clus_metric, B = n_B, spaceH0 = ref_dist,
                             d_space = d,
                             FUN = FUN, center_F_marg = TRUE, max_iter = max_iter, tol = tol,
                             foldid = foldid)
  list(gap_stat = gap_stat, er_init = er_init)
  
}


clusGap_mixpfc <- function (x, y, K.max, gams_init_list,
                            lam_fact_list, within_clus_metric = "ss",
                            B = 100, d.power = 1,
                            spaceH0 = c("scaledPCA", "original"),
                            d_space = NULL,
                            verbose = interactive(), ...) {
  # input:
  # within_clus_metric: "ss", sum of square, "angle", sum of angles
  #                     "ss_xy", sum of square on matrix xy
  stopifnot(length(dim(x)) == 2, K.max >= 2,
            (n <- nrow(x)) >= 1, ncol(x) >= 1)
  if(B != (B. <- as.integer(B)) || (B <- B.) < 0)
    stop("'B' has to be a non-negative integer")
  cl. <- match.call()
  
  if(is.data.frame(x))
    x <- as.matrix(x)
  ii <- seq_len(n)
  W.k <- function(X, Y, kk) {
    clus_space <- if(kk > 1) {
      mixpfc_fit <- mixPFC_K(X, Y, 
                             n_components = kk, gams_init = gams_init_list[[kk-1]],
                             weight_init = colMeans(gams_init_list[[kk-1]]),
                             lam_factor = lam_fact_list[[kk-1]],
                             ...)
      clus <- est_component(mixpfc_fit$gams)
      space <- mixpfc_fit$B_mats
      row_ind_hat <- lapply(space, get_ind_row)
      row_ind_hat_union <- Reduce(union, row_ind_hat) # the important variable set
      if (length(row_ind_hat_union) < ncol(X)) {
        # refit with important variables
        refit <- mixPFC_K(X[, row_ind_hat_union, drop = FALSE], Y, 
                          n_components = kk, gams_init = mixpfc_fit$gams,
                          weight_init = colMeans(mixpfc_fit$gams),
                          lam_factor = rep(0, kk),
                          ...)
        clus <- est_component(refit$gams)
        space <- refit$B_mats
      } 
      list(clus = clus, space = space, act_var = row_ind_hat_union)
    } else {
      clus <- rep.int(1L, nrow(X))
      space <- diag(ncol(X))
      row_ind_hat_union <- seq_len(ncol(X))
      list(clus = clus, space = space, act_var = row_ind_hat_union)
    }
    
    clus <- clus_space$clus
    act_var <- clus_space$act_var
    # cat(act_var, "\n")
    X <- X[, act_var, drop = FALSE] # redefine X to only include the important variables
    
    if (within_clus_metric == "ss") {
      val <- 0.5* sum(vapply(split(ii, clus),
                             function(I) { xs <- X[I,, drop=FALSE]
                             sum(dist(xs)^d.power/nrow(xs)) }, 0.))
    } else if (within_clus_metric == "ss_xy") {
      XY <- X * matrix(Y, nrow = nrow(X), ncol = ncol(X))
      val <-  0.5* sum(vapply(split(ii, clus),
                              function(I) { xs <- XY[I,, drop=FALSE]
                              sum(dist(xs)^d.power/nrow(xs)) }, 0.))
    } else if (within_clus_metric == "angle") {
      val <-  0.5* sum(vapply(split(ii, clus),
                              function(I) { 
                                xs <- X[I,, drop=FALSE]
                                angle_mat <- pairwise_angle_matrix(xs) * pi / 180
                                sum(angle_mat^d.power/nrow(xs)) }, 
                              0.))
    } else if (within_clus_metric == "ss_proj_orth") {
      # if we project X onto the subspace that is orthogonal to beta
      if (is.null(d_space)) {
        stop("clusGap_mixpfc: d_space is null!")
      }
      
      if (kk == 1) {
        X_trasformed <- X
      } else {
        # assume all space have the same dimension
        X_trasformed <- matrix(0, nrow = nrow(X), ncol = ncol(X)-d_space)
        B_hat <- clus_space$space
        beta0_hat <- vector("list", length = kk)
        for (j in seq_len(kk)) {
          beta_hat <- svd(B_hat[[j]])$u[, 1:d_space, drop = FALSE] # assume same dimension
          beta0_hat[[j]] <- orthogonal_complement(beta_hat)
          clus_j <- which(clus == j)
          X_trasformed[clus_j, ] <- X[clus_j, ,drop=FALSE] %*% beta0_hat[[j]]
        }
      }
      
      val <- 0.5* sum(vapply(split(ii, clus),
                             function(I) { xs <- X_trasformed[I,, drop=FALSE]
                             sum(dist(xs)^d.power/nrow(xs)) }, 0.))
    }
    val
  }
  logW <- E.logW <- SE.sim <- numeric(K.max)
  if(verbose) cat("Clustering k = 1,2,..., K.max (= ",K.max,"): .. ", sep='')
  for(k in 1:K.max)
    logW[k] <- log(W.k(x, y, k))
  if(verbose) cat("done\n")
  
  if (B == 0) {
    # B = 0 means only compute W
    return(logW)
  }
  
  spaceH0 <- match.arg(spaceH0)
  ## Scale 'x' into hypercube -- later fill with H0-generated data
  xs <- scale(x, center=TRUE, scale=FALSE)
  m.x <- rep(attr(xs,"scaled:center"), each = n) # for back-trafo later
  switch(spaceH0,
         "scaledPCA" =
           {
             ## (These & (xs,m.x) above basically do stats:::prcomp.default()
             V.sx <- svd(xs, nu=0)$v
             xs <- xs %*% V.sx # = transformed(x)
           },
         "original" = {}, # (do nothing, use 'xs')
         ## otherwise
         stop("invalid 'spaceH0':", spaceH0))
  
  rng.x1 <- apply(xs, 2L, range)
  rng.y1 <- range(y)
  logWks <- matrix(0, B, K.max)
  if(verbose) cat("Bootstrapping, b = 1,2,..., B (= ", B,
                  ")  [one \".\" per sample]:\n", sep="")
  for (b in 1:B) {
    ## Generate "H0"-data as "parametric bootstrap sample" :
    z1 <- apply(rng.x1, 2,
                function(M, nn) runif(nn, min=M[1], max=M[2]),
                nn=n)
    z <- switch(spaceH0,
                "scaledPCA" = tcrossprod(z1, V.sx), # back transformed
                "original" = z1
    ) + m.x
    ## Generate "H0"-response
    y1 <- runif(n, min = rng.y1[1], max = rng.y1[2])
    for(k in 1:K.max) {
      logWks[b,k] <- log(W.k(z, y1, k))
    }
    if(verbose) cat(".", if(b %% 50 == 0) paste(b,"\n"))
  }
  if(verbose && (B %% 50 != 0)) cat("",B,"\n")
  E.logW <- colMeans(logWks)
  SE.sim <- sqrt((1 + 1/B) * apply(logWks, 2, var))
  structure(class = "clusGap",
            list(Tab = cbind(logW, E.logW, gap = E.logW - logW, SE.sim),
                 ## K.max == nrow(T)
                 call = cl., spaceH0=spaceH0,
                 n = n, B = B))
}




# Function to compute pairwise angles between rows of a matrix
pairwise_angle_matrix <- function(X) {
  n <- nrow(X)
  angle_matrix <- matrix(NA, n, n)
  for (i in 1:n) {
    for (j in i:n) {
      if (i == j) {
        angle_matrix[i, j] <- 0
        next
      }
      angle <- compute_angle(X[i, ], X[j, ])
      # angle <- ifelse(angle > 90, 180 - angle, angle)
      angle_matrix[i, j] <- 0
      angle_matrix[j, i] <- angle  # Since the matrix is symmetric
    }
  }
  angle_matrix
}






















