# Construct AR matrix
AR <- function(rho, p){
  m <- matrix(0, p, p)
  for (i in 1:p){
    for (j in 1:p){
      m[i,j] <- rho**(abs(i-j))
    }
  }
  return(m)
}

# Cut small values in a matrix to zero. 
cut_mat <- function(Beta, thrd, rank){
  l <- length(Beta)
  for (i in 1:l){
    if(is.null(Beta[[i]])) next
    mat <- as.matrix(Beta[[i]])
    nobs <- nrow(mat)
    nvars <- ncol(mat)
    r <- rank[i]
    if(r == 0){
      Beta[[i]] <- matrix(0, nobs, nvars)
    }else{
      vec <- as.vector(mat)
      vec[abs(vec) < thrd] <- 0
      Beta[[i]] <- matrix(vec, nobs, nvars)
    }
  }
  return(Beta)
}

# Evaluation based on distance correlation. (Szekely et al., 2007)
eval_dc <- function(Beta, x, y, type, gam = NULL){
  if(!is.list(Beta)){Beta <- list(Beta)}
  l <- length(Beta)
  result <- sapply(seq_len(l), function(i){
    if(is.null(Beta[[i]])){
      NA
    }else{
      mat <- as.matrix(Beta[[i]])
      if (type != "mixpfc") {
        dcor(x %*% mat, y)
      } else {
        dcor(diag(sqrt(gam)) %*% x %*% mat, diag(sqrt(gam)) %*% y)
      }
    }
  })
  return(result)
}

# Compute M and U matrices from observation data.
#########
# Input:
# x: n x p observation matrix for predictor.
# y: n-dimensional observation vector for response.
# yclass: Discretized response taking values in 1,...,H.
# type: Specifying the specific SEAS method. "sir" means SEAS-SIR, "intra" means SEAS-Intra and "pfc" means SEAS-PFC.
# FUN: the user-specified function f in SEAS-PFC. The default is f(y) = (y, y^2, y^3).
# categorical: A logical value indicating whether y is categorical.
# H: The number of slices. The default value is 5.
# gams: a length n vector, the posterior probabilities of each observations for 
#       a component
# sig_fit: boolean, if return the matrix sigma_fit
# center_F_marg: boolean, if we center the F matrix marginally intead of 
#               center F using gam
########################## following paras are for mixpfc ######################
# gam: n dimensional vector, which is the j-th column of gams (n*n_components matrix)

MU <- function(x, y, yclass=NULL, type='sir', FUN = NULL, categorical = FALSE, H = 5,
               gam = NULL, sig_fit = FALSE, center_F_marg = FALSE){
  if(is.null(yclass)){ # Construct the discretized response
    if(categorical == FALSE){
      ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
      yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
      nclass <- as.integer(length(unique(yclass)))
    }
    else if(categorical == TRUE){
      yclass <- y
    }
  }
  cls <- sort(unique(yclass))
  nclass <- length(cls)
  nobs <- as.integer(dim(x)[1])
  nvars <- as.integer(dim(x)[2])
  prior <- sapply(cls, function(i){mean(yclass == i)})
  mu <- colMeans(x)
  x_c <- x - matrix(mu, nobs, nvars, byrow = TRUE) # centered predictor
  M <- crossprod(x_c/sqrt(nobs)) # sample covariance of X
  
  if(type == 'sir'){
    U <- matrix(0, nvars, nclass)
    for (i in 1:nclass){
      U[, i] <- colMeans(x_c[yclass == cls[i],, drop=FALSE])
    }
    sigma_fit <- Fmat_c <- NA
  }else if(type == 'intra'){
    y_c <- y - mean(y)
    U <- matrix(0, nvars, nclass)
    lb <- quantile(y_c, 0.1)[[1]]
    ub <- quantile(y_c, 0.9)[[1]]
    y_c <- sapply(y_c, cut_func, lb = lb, ub = ub)
    for (i in 1:nclass){
      y_copy <- y_c
      y_copy[yclass!=cls[i]] <- 0
      U[, i] <- (1/nobs) * t(x_c) %*% (y_copy - mean(y_copy))
    }
  }else if(type == 'pfc'){
    if(is.null(FUN)) Fmat <- cbind(y, y^2, y^3) # the default function
    else Fmat <- t(sapply(y, FUN))
    if (nrow(Fmat) == 1) Fmat <- t(Fmat) # if the input function is one function
    Fmat_mean <- colMeans(Fmat)
    Fmat_c <- Fmat - matrix(Fmat_mean, nrow(Fmat), ncol(Fmat), byrow = TRUE) # centered function f
    lb <- apply(Fmat_c, 2, quantile, 0.1)
    ub <- apply(Fmat_c, 2, quantile, 0.9)
    for(i in 1:NCOL(Fmat_c)){
      Fmat_c[, i] <- sapply(Fmat_c[,i], cut_func, lb[i], ub[i])
    }
    U <- (1/nobs)*(t(x_c) %*% Fmat_c)
    
    # generate sigma fit
    if (sig_fit) {
      P_F <- Fmat_c %*% ginv(t(Fmat_c) %*% Fmat_c) %*% t(Fmat_c)
      sigma_fit <- (1/nobs) * t(x_c) %*% P_F %*% x_c
    } else {
      sigma_fit <- NA
    }
  } else if (type == 'mixpfc') {
    #Fmat_c <- compute_Fmat(y, FUN)
    if(is.null(FUN)) Fmat <- cbind(y, y^2, y^3) # the default function
    else Fmat <- t(sapply(y, FUN))
    if (nrow(Fmat) == 1) Fmat <- t(Fmat)
    

    ### compute Fmat_mean using gam
    Fmat_mean <- rep(0, ncol(Fmat))
    
    if (center_F_marg) { ## if center_F_marg is true
      Fmat_mean <- colMeans(Fmat)
    } else {            ## otherwise use gam to center F
      for (i in seq_len(nrow(Fmat))) {
        Fmat_mean <- Fmat_mean + gam[i] * Fmat[i, ]
      }
      Fmat_mean <- Fmat_mean / sum(gam)
    }
    ###
    Fmat_c <- Fmat - matrix(Fmat_mean, nrow(Fmat), ncol(Fmat), byrow = TRUE) # centered function f
    lb <- apply(Fmat_c, 2, quantile, 0.1)
    ub <- apply(Fmat_c, 2, quantile, 0.9)
    for(i in 1:NCOL(Fmat_c)){
      Fmat_c[,i] <- sapply(Fmat_c[,i], cut_func, lb[i], ub[i])
    }
    mu <- rep(0, nvars)
    for (i in seq_len(nobs)) {
      mu <- mu + gam[i] * x[i, ]
    }
    mu <- mu / sum(gam)
    x_c <- x - matrix(mu, nobs, nvars, byrow = TRUE) # centered predictor
    M <- (1/sum(gam)) * (t(x_c) %*% diag(gam) %*% x_c)
    U <- (1/sum(gam)) * (t(x_c) %*% diag(gam) %*% Fmat_c)
    if (sig_fit) {
      #Fmat_c_gam <- diag(sqrt(gam)) %*% Fmat_c
      #P_F <- Fmat_c %*% ginv(t(Fmat_c) %*% Fmat_c) %*% t(Fmat_c)
      #cat("Num of NA in gam: ", sum(is.na(gam)))
      #print(t(Fmat_c) %*% diag(gam) %*% Fmat_c)
      P_F <- Fmat_c %*% ginv(t(Fmat_c) %*% diag(gam) %*% Fmat_c) %*% t(Fmat_c)
      sigma_fit <- (1/sum(gam)) * (t(x_c) %*% diag(gam) %*% P_F %*% diag(gam) %*% x_c)
    } else {
      sigma_fit <- NA
    }
  }
  list(M = M, U = U, nclass = nclass, prior=prior, sigma_fit = sigma_fit, 
       Fmat_c = Fmat_c, x_c = x_c)
}

# Cut extreme values in the samples. This function is used in MU function.
cut_func <- function(x, lb, ub){
  if(x < lb){
    return(lb)
  } else if(x > ub){
    return(ub)
  } else{
    return(x)
  }
}

# Estimate the rank of a matrix.
rank_func <- function(B, thrd){
  d <- svd(B)$d
  r <- sum(d >= thrd)
  return(r)
}

# Subspace distance, defined in (19)
subspace <- function(A,B){
  if(is.vector(A)) A <- as.matrix(A)
  if(is.vector(B)) A <- as.matrix(B)
  Pa <- qr.Q(qr(A))
  Pa <- Pa %*% t(Pa)
  Pb <- qr.Q(qr(B))
  Pb <- Pb %*% t(Pb)
  d <- dim(A)[2]
  return(norm(Pa-Pb, type="F")/sqrt(2*d))
}

## ------------------------------------------------ ##
## The utility functions imported from R package 'msda'.
## These functions are used in 'msda' and 'cv.msda' functions.
formatoutput <- function(fit, maxit, pmax, p, H) {
  nalam <- fit$nalam
  ntheta <- fit$ntheta[seq(nalam)]
  nthetamax <- max(ntheta)
  lam <- fit$alam[seq(nalam)]
  theta_vec <- fit$theta
  errmsg <- err(fit$jerr, maxit, pmax)  ### error messages from fortran
  switch(paste(errmsg$n), `1` = stop(errmsg$msg, call. = FALSE), `-1` = cat(errmsg$msg))
  if(nthetamax > 0){
    ja <- fit$itheta[seq(nthetamax)]
    theta <- lapply(seq_len(nalam), function(i){
      tmp <- theta_vec[(pmax * H * (i-1) + 1):(pmax * H * i)]
      a <- matrix(tmp, pmax, H, byrow = TRUE)[seq(nthetamax), , drop = FALSE]
      theta_i <- matrix(0, p, H)
      theta_i[ja,] <- a
      theta_i
    })
  }
  else{
    theta <- lapply(seq(nalam), function(x){matrix(0, p, H)})
  }
  list(theta = theta, lambda = lam)
}

err <- function(n, maxit, pmax) {
  if (n == 0) 
    msg <- ""
  if (n > 0) {
    # fatal error
    if (n < 7777) 
      msg <- "Memory allocation error; contact package maintainer"
    if (n == 10000) 
      msg <- "All penalty factors are <= 0"
    n <- 1
    msg <- paste("in the fortran code -", msg)
  }
  if (n < 0) {
    # non fatal error
    if (n > -10000) 
      msg <- paste("Convergence for ", -n, "th lambda value not reached after maxit=", maxit, " iterations; solutions for larger lambdas returned.\n", sep = "")
    if (n < -10000) 
      msg <- paste("Number of nonzero coefficients along the path exceeds pmax=", pmax, " at ", -n - 10000, "th lambda value; solutions for larger lambdas returned.\n", sep = "")
    if (n < -20000) 
      msg <- paste("Number of nonzero coefficients along the path exceeds dfmax=", pmax, " at ", -n - 20000, "th lambda value; solutions for larger lambdas returned.\n", sep = "")
    n <- -1
  }
  list(n = n, msg = msg)
}

lamfix <- function(lam){
  llam <- log(lam)
  if(length(llam) >= 3){lam[1] <- exp(2 * llam[2] - llam[3])}
  lam
}

################################################################################
################################################################################

################## utilities functions for mixPFC ##############################

################################################################################
################################################################################

# Compute the sample covariance matrix
# Example:
# n <- 5  # Number of samples
# p <- 3  # Number of variables
# set.seed(2023)
# X <- matrix(rnorm(n * p), nrow = n, ncol = p)
# cov_mat(X)
cov_mat <- function(X) {
  n <- as.integer(nrow(X))
  p <- as.integer(ncol(X))
  mu <- colMeans(X)
  X_c <- X - matrix(mu, n, p, byrow = TRUE) # centered
  Sigma <- crossprod(X_c) / n
  Sigma
}

cov_fit <- function(X, Fmat) {
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Fmat)
  mu <- colMeans(X)
  X_c <- X - matrix(mu, n, p, byrow = TRUE) # centered
  
  F_mean <- colMeans(Fmat)
  Fmat_c <- Fmat - matrix(F_mean, n, q, byrow = TRUE) # centered
  P_F <- Fmat_c %*% ginv(t(Fmat_c) %*% Fmat_c) %*% t(Fmat_c)
  Sigma_fit <- t(X_c) %*% P_F %*% X_c / n
  Sigma_fit
}

center_mat <- function(X) {
  # center matrix X columnwise
  n <- nrow(X)
  p <- ncol(X)
  mu <- colMeans(X)
  X_c <- X - matrix(mu, n, p, byrow = TRUE) # centered
  X_c
}

is_psd <- function(X) {
  # true is matrix X is psd, false elsewise
  if (min(svd(X)$d) < 1e-5) return(FALSE)
  else return(TRUE)
}





# generate initial gams with only 2 cluster
gen_gams_init <- function(component_ind, true_percent = 0.8, true_comp_ind = TRUE) {
  # component_ind: length n vector, the true component id
  # true_percent: a real number between 0 and 1, used only if true_comp_ind is TRUE
  # true_comp_ind: TRUE if component_ind is the true ind of data, in this case
  #                n*true_percent samples would have the true id
  #                If true_comp_ind is FALSE, we will use component_ind to 
  #                generate gams_init
  n <- length(component_ind)
  n_components <- length(unique(component_ind))
  if (n_components != 2) cat("number of components must be 2!, n_components=", n_components)
  n_components <- 2
  gams_init <- matrix(0, n, n_components)
  if (true_comp_ind) { # if component_ind is the true component id
    n_true <- floor(n * true_percent)
    for (i in seq_len(n)) {
      if (i <= n_true) {
        gams_init[i, component_ind[i]] <- 1
      }
      else {
        if (component_ind[i] == 1)
          gams_init[i, 2] <- 1
        else
          gams_init[i, 1] <- 1
      }
    }
  } else { # if component_ind is not the true component id
    for (i in seq_len(n)) {
      gams_init[i, component_ind[i]] <- 1
    }
  }
  
  gams_init
}

# generate initial gams with more than 2 clusters
gen_gams_init_K <- function(component_ind, true_percent = 0.7, true_comp_ind = TRUE) {
  # component_ind: length n vector, the true component id
  # true_percent: a real number between 0 and 1, used only if true_comp_ind is TRUE
  # true_comp_ind: TRUE if component_ind is the true ind of data, in this case
  #                n*true_percent samples would have the true id
  #                If true_comp_ind is FALSE, we will use component_ind to 
  #                generate gams_init
  
  # DO NOT SET SEED FOR THIS FUNCTION!!!
  # typically, we use this function after the data generating function
  # if we set seed here, then the data gen function will generate SAME data 
  # for different repetition
  # set.seed(2023) 
  n <- length(component_ind)
  n_components <- length(unique(component_ind))
  gams_init <- matrix(0, n, n_components)
  if (true_comp_ind) { # if component_ind is the true component id
    n_true <- floor(n * true_percent)
    # randomly select 
    true_indices <- sample(n, n_true)
    
    for (i in true_indices) {
      gams_init[i, component_ind[i]] <- 1
    }
    
    remaining_indices <- setdiff(1:n, true_indices)
    for (i in remaining_indices) {
      # randomly select a column index between 1 and K
      
      random_cluster <- sample(setdiff(1:n_components, component_ind[i]), 1)
      gams_init[i, random_cluster] <- 1
    }
    
  } else { # if component_ind is not the true component id
    for (i in seq_len(n)) {
      gams_init[i, component_ind[i]] <- 1
    }
  }
  
  gams_init
}

# compute the posterior probability of data (X, Y) for general covariance matrix
# we have another function for isotonic mixPFC with eta in "init_utility.R"
compute_post_without_gam <- function(X, Y, B1, B2, weight_pi, A, Delta_inv,
                                     FUN = NULL, mu = NULL) {
  # Input:
  # X, Y: data
  # B1, B2: p*q matrices
  # A: a list with length = 2, each element is a 2q*q matrix
  # Delta_inv: 2q*2q matrix
  # weight_pi: weights for each cluster
  # FUN: function to compute Fmat
  # mu: a list with length = 2, cluster mean
  
  n <- nrow(X) # sample size
  n_components <- 2
  gams <- matrix(nrow = n, ncol = n_components)
  #B <- cbind(B1, B2) # p*2q matrix
  
  X_c <- vector("list", n_components) # a list, centered by mu
  for (j in seq_len(n_components)) {
    X_c[[j]] <- X - matrix(mu[[j]], nrow(X), ncol(X), byrow = TRUE)
  }
  
  Z <- cbind(X_c[[1]] %*% B1, X_c[[2]] %*% B2) # n*2q matrix
  Fmat <- compute_Fmat(Y, FUN = FUN)
  Fmat_c <- center_Fmat(Fmat)
  for (i in seq_len(n)) {
    num_exp <- rep(0, n_components) # record the number in the exp (without the nagtive sign)
    for (j in seq_len(n_components)) {
      temp <- matrix(Z[i,] - matrix(Fmat_c[i, ], nrow = 1)%*%t(A[[j]]), nrow = 1)
      #gams[i, j] <- weight_pi[j] *
      #  exp(-0.5*temp %*% Delta_inv %*% t(temp))
      num_exp[j] <- 0.5*temp %*% Delta_inv %*% t(temp)
    }
    gams[i, 1] <- weight_pi[1] / (weight_pi[1] + weight_pi[2]*exp(num_exp[1]-num_exp[2]))
    gams[i, 2] <- weight_pi[2] / (weight_pi[1]*exp(num_exp[2]-num_exp[1]) + weight_pi[2])
  }
  gams
}



# visualize data
visul_data <- function(X, Y, beta, ind, d, sample_fraction = 1) {
  # this function first project X onto subspaces spanned by each beta
  # and the draw a figure with different colors
  # X: n*p dimension matrix
  # Y: n*1 vector
  # beta: a list, each is a basis of the central subspace of dimension p*d
  #       assuming d = 2
  # ind: a length-n vector consisting of indicator of which components a sample
  #      belongs
  # sample_fraction: a number between 0 and 1, control how many sample we use
  
  # Define RGB colors for each level of the "Group" factor
  
  n <- nrow(X) # sample size
  num_keep <- round(n * sample_fraction)
  sample_ind <- sample(seq_len(n), size = num_keep, replace = FALSE)
  X <- X[sample_ind, ]
  Y <- Y[sample_ind]
  ind <- ind[sample_ind]

  #color1 <- "#FF0000" # red 
  #color1 <- "#00008B" # dark blue
  color1 <- "#d7191c"
  color2 <- "#2b83ba"
  colors <- c(color1, color2)
  shapes <- c(16, 17) # shapes of the points, 16 is circle and 17 is triangle
  if (d == 1) {
    plot_title1 <- expression(paste("Sufficient Summary Plot", 
                                " (", "d=1", ")", sep = ""))
    plot_title2 <- plot_title3 <- plot_title1
    #dirc1 <- dirc2 <- expression(paste("Predictor Projected onto ", 
    #                                   S[1], " (=", S[2], ")", sep = ""))
    dirc1 <- dirc2 <- latex2exp::TeX(r'(Predictor Projected onto $S_1$ ($=S_2$))')
    
  } else if (d == 2) {
    plot_title1 <- expression(paste("Sufficient Summary for Mixture 1", 
                            " (", d[1], "=", 1, ")", sep = ""))
    plot_title2 <- expression(paste("Sufficient Summary for Mixture 2", 
                                    " (", d[2], "=", 1, ")", sep = ""))
    plot_title3 <- expression(paste("Predictor Reduction", 
                                    " (", "d=2", ")", sep = ""))
    # dirc1 <- expression(paste("Predictor Projected onto ", 
    #                           S[1], sep = ""))
    # dirc2 <- expression(paste("Predictor Projected onto ", 
    #                           S[2], sep = ""))
    dirc1 <- latex2exp::TeX(r'(Predictor Projected onto $S_1$ ($\perp S_2$))')
    dirc2 <- latex2exp::TeX(r'(Predictor Projected onto $S_2$ ($\perp S_1$))')
  }
  
  # Project data onto the two subspaces
  X_proj <- cbind(X %*% beta[[1]], X %*% beta[[2]])
  
  # Create a data frame
  df <- data.frame(Direction1 = X_proj[, 1], 
                   Direction2 = X_proj[, 2],
                   Y = Y,
                   Mixture = factor(ind))
  
  
  
  # Create a scatter plot with different colors for each matrix
  p1 <- ggplot(df, aes(x = Direction1, y = Y, color = Mixture,
                       shape = Mixture)) + 
    geom_point(size=4.5, stroke = 1.5) + # size = 6, stroke = 3; size = 3.3, stroke = 1.5
    scale_color_manual(values = colors) +
    scale_shape_manual(values = c(21, 17)) +
    #scale_size_manual(values = c(4, 4)) +
    labs(x = dirc1, y = "Response") +
    #ggtitle(plot_title1) + 
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_blank(),
          axis.title.x = element_text(size = 26),
          axis.title.y = element_text(size = 26),
          legend.position = "none")
  
  p2 <- ggplot(df, aes(x = Direction2, y = Y, color = Mixture,
                       shape = Mixture)) + 
    geom_point(size=4.5, stroke = 1.5) + # size = 6, stroke = 3; size = 3.3, stroke = 1.5
    scale_color_manual(values = colors) +
    scale_shape_manual(values = c(1, 17)) +
    labs(x = dirc2, y = "Response") +
    #ggtitle(plot_title2) + 
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_blank(),
          axis.title.x = element_text(size = 26),
          axis.title.y = element_text(size = 26),
          legend.position = "none")
  
  p3 <- ggplot(df, aes(x = Direction1, y = Direction2, color = Mixture,
                       shape = Mixture)) + 
    geom_point(size=1.5) +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = c(16, 4)) +
    labs(x = dirc1, y = dirc2) +
    #ggtitle(plot_title3) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_blank(),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 18))
  
  
  list(p1 = p1, p2 = p2, p3 = p3)
}



visul_data_3D_2 <- function(X, Y, beta, ind) {
  # this function first project X onto subspaces spanned by each beta
  # and then draws a 3D scatter plot with different colors
  # X: n*p dimension matrix
  # Y: length n vector
  # beta: a list, each is a basis of the central subspace of dimension p*d
  #       assuming d = 1
  # ind: a vector specifying the group membership for each observation
  
  # Project data onto the two subspaces
  X_proj <- cbind(X %*% beta[[1]], X %*% beta[[2]])
  
  # point_shape <- rep("circle", length(ind))
  # point_shape[ind==2] <- "x"
  
  # Create a data frame
  df <- data.frame(Direction1 = X_proj[, 1], 
                   Direction2 = X_proj[, 2],
                   Y = Y,
                   Group = factor(ind))
  
  #Create the 3D scatter plot using plotly
  p <- plot_ly(data = df, x = ~Direction1, y = ~Direction2, z = ~Y, color = ~Group,
          colors = c("#66C2A5", "#FC8D62"),
          symbol = ~Group,
          symbols = c("circle", "x"),
          type = "scatter3d", mode = "markers",
          # symbol = point_shape,
          marker = list(size = 4)) %>%
    layout(scene = list(xaxis = list(title = "Predictor Reduction by S<sub>1</sub>",
                                     titlefont = list(size = 19),  # Set x-axis label font size
                                     showticklabels = FALSE),
                        yaxis = list(title = "Predictor Reduction by S<sub>2</sub>",
                                     titlefont = list(size = 19), 
                                     showticklabels = FALSE),
                        zaxis = list(title = "Response Y ",
                                     titlefont = list(size = 19), 
                                     showticklabels = FALSE),
                        aspectmode = "cube"),
           title = "3D Sufficient Summary Plot (d=2, d<sub>1</sub>=d<sub>2</sub>=1)",
           legend = list(title = list(text="Mixture"), traceorder = "normal",
                         bgcolor = "#FFFFFF", bordercolor = "#FFFFFF",
                         itemsizing = "constant",
                         font = list(size = 19),  # Set legend text size
                         x=0.85, y=0.65))
  

  # p <- plotly_build(p)
  # p$x$data[[1]]$marker$symbol <- "x"
  p
}



center_X_F_by_gam <- function(X, Fmat, gam, center_F_marg = FALSE) {
  # X:n * p
  # F: n* q
  # gam: length n vector
  # center_F_marg: boolean, if we center the F matrix marginally intead of 
  #               center F using gam
  nobs <- nrow(X)
  nvars <- ncol(X)
  mu <- rep(0, nvars)
  for (i in seq_len(nobs)) {
    mu <- mu + gam[i] * X[i, ]
  }
  mu <- mu / sum(gam)
  X_c <- X - matrix(mu, nobs, nvars, byrow = TRUE) # centered predictor
  
  Fmat_mean <- rep(0, ncol(Fmat))
  if (center_F_marg) { # if center F marginally
    Fmat_mean <- colMeans(Fmat)
  } else {             # if center F using gam
    for (i in seq_len(nrow(Fmat))) {
      Fmat_mean <- Fmat_mean + gam[i] * Fmat[i, ]
    }
    Fmat_mean <- Fmat_mean / sum(gam)
  }

  ###
  Fmat_c <- Fmat - matrix(Fmat_mean, nrow(Fmat), ncol(Fmat), byrow = TRUE) # centered function f
  lb <- apply(Fmat_c, 2, quantile, 0.1)
  ub <- apply(Fmat_c, 2, quantile, 0.9)
  for(i in 1:NCOL(Fmat_c)){
    Fmat_c[,i] <- sapply(Fmat_c[,i], cut_func, lb[i], ub[i])
  }
  list(X_c = X_c, Fmat_c = Fmat_c)
}


center_X_by_gam <- function(X, gam) {
  # X:n * p
  # gam: length n vector

  nobs <- nrow(X)
  nvars <- ncol(X)
  mu <- rep(0, nvars)
  for (i in seq_len(nobs)) {
    mu <- mu + gam[i] * X[i, ]
  }
  mu <- mu / sum(gam)
  X_c <- X - matrix(mu, nobs, nvars, byrow = TRUE) # centered predictor
  
  X_c
}

# get active predictor, the indices of zero rows of estimated beta
get_ind_row <- function(beta_hat) {
  which(apply(beta_hat, 1, function(x) any(x != 0)))
}

# compute true positive rate (TPR) and false positive rate (FPR)
# let A be the true active set, TPR = |A intersection A^hat| / |A|
# and FPR = |A_hat intersection A^c| / |A^c|
compute_TPR <- function(A_hat, A) {
  # A_hat: a vector containing estimated ind of active variables
  # A: a vector containing true ind of active variables
  sum(A_hat %in% A) / length(A)
}

compute_FPR <- function(A_hat, A, p) {
  # A_hat: a vector containing estimated ind of active variables
  # A: a vector containing true ind of active variables
  # p: dimension of predictors
  ifelse(p == length(A), 0, sum(!(A_hat %in% A)) / (p - length(A)))
}


est_component <- function(gams) {
  #gams: n_sample * n_component
  n <- nrow(gams)
  est_comp <- rep(0, n)
  for (i in seq_len(n)) {
    est_comp[i] <- which.max(gams[i, ])[1]
  }
  
  est_comp
}

get_est_sub_from_B <- function(B_hat, d1, d2) {
  # B_hat: a list of estimated matrices, after soft-threshold of the singular
  #       values we can obtain the estimated subspaces
  # d1, d2: true dimension of the subspaces
  #         if d1 = d2, the order of beta_hat matters 
  
  # consider the combination of 2 clusters and 2 dimension
  # first group
  beta_hat1 <- svd(B_hat[[1]])$u[, 1:d1, drop = FALSE]
  beta_hat2 <- svd(B_hat[[2]])$u[, 1:d2, drop = FALSE]
  # second group
  beta_hat3 <- svd(B_hat[[2]])$u[, 1:d1, drop = FALSE]
  beta_hat4 <- svd(B_hat[[1]])$u[, 1:d2, drop = FALSE]
  
  list(list(beta_hat1, beta_hat2, id=c(1, 2)),
       list(beta_hat3, beta_hat4, id=c(2, 1)))
}

sub_dist_error_helper <- function(beta, beta_hat) {
  # we assume the order of beta_hat is correct
  # thus just return the two distance
  # input:
  # beta: a list with 2 elements, the true subspaces
  # beta_hat: a list with 2 elements, the estimated matrices whose column
  #           spaces estimate the true subspaces
  # output:
  #   id: a length 2 vector, beta_hat[[id[i]]] is the estimate of beta[[i]]
  sub_dist11 <- subspace(beta[[1]], beta_hat[[1]])
  #sub_dist12 <- subspace(beta[[1]], beta_hat[[2]])
  
  #sub_dist21 <- subspace(beta[[2]], beta_hat[[1]])
  sub_dist22 <- subspace(beta[[2]], beta_hat[[2]])
  
  # if ((sub_dist11 + sub_dist22) <= (sub_dist12 + sub_dist21)) {
  #   dist1_error <- sub_dist11
  #   dist2_error <- sub_dist22
  #   id <- c(1, 2) # id: (1, 2) means beta_hat[[1]] is the estimator of beta[[1]]
  #   beta_est <- beta_hat
  # } else {
  #   dist1_error <- sub_dist12
  #   dist2_error <- sub_dist21
  #   id <- c(2, 1) # id: (2, 1) means beta_hat[[2]] is the estimator of beta[[1]]
  #   beta_est <- list(beta_hat[[2]], beta_hat[[1]])
  # }
  # list(dist_error = c(dist1_error, dist2_error),
  #      id = id, beta_est = beta_est)
  c(sub_dist11, sub_dist22)
}

sub_dist_error <- function(beta, B_hat) {
  # beta: a list with 2 elements, the true subspaces
  # B_hat: a list of estimated matrices, after soft-threshold of the singular
  #       values we can obtain the estimated subspaces
  
  d1 <- ncol(beta[[1]]) # dimension of first subspace
  d2 <- ncol(beta[[2]]) # dimension of second subspace
  
  beta_hat <- get_est_sub_from_B(B_hat, d1, d2)
  beta_hat1 <- list(beta_hat[[1]][[1]], beta_hat[[1]][[2]])
  beta_hat2 <- list(beta_hat[[2]][[1]], beta_hat[[2]][[2]])
  
  dist1 <- sub_dist_error_helper(beta, beta_hat1)
  dist2 <- sub_dist_error_helper(beta, beta_hat2)
  
  if (sum(dist1) <= sum(dist2)) {
    return(list(dist_error=dist1, id=beta_hat[[1]]$id))
  } else {
    return(list(dist_error=dist2, id=beta_hat[[2]]$id))
  }
  
}

sub_dist_error_K <- function(beta, beta_hat) {
  # different than function `sub_dist_error` for two clusters
  # the second input parameter is beta_hat, which means the matrix B_hat
  # is already cut into beta_hat
  # Input: 
  # beta: a list with K elements, the basis of true subspaces
  # beta_hat: a list with K element, the estimated basis of true subspaces
  n_components <- length(beta)
  dist_error <- rep(Inf, n_components)
  per_K <- combinat::permn(n_components) # permuation of c(1:K)
  
  for (i in seq_along(per_K)) {
    # beta_hat for the i-th permutation
    beta_hat_curr <- beta_hat[per_K[[i]]]
    # then compute subspace distance
    dist_error_curr <- rep(NA, n_components)
    for (k in seq_len(n_components)) {
      dist_error_curr[k] <- subspace(beta[[k]], beta_hat_curr[[k]])
    }
    
    if (sum(dist_error_curr) < sum(dist_error)) {
      dist_error <- dist_error_curr
      matched_perm <- per_K[[i]] # the matched permutation
      beta_hat_matched <- beta_hat_curr 
    }
  }
  
  list(dist_error = dist_error, matched_perm = matched_perm,
       beta_hat_matched = beta_hat_matched)
  
}

get_dcor <- function(x, y, comp_hat, beta_hat) {
  # x, y: data
  # comp_hat: estimated component
  # beta_hat: a list of estimated matrices
  
  ind1 <- comp_hat == 1
  dcor11 <- dcor(x[ind1, ] %*% beta_hat[[1]], y[ind1]) 
  dcor22 <- dcor(x[!ind1, ] %*% beta_hat[[2]], y[!ind1])
  
  dcor12 <- dcor(x[ind1, ] %*% beta_hat[[2]], y[ind1]) 
  dcor21 <- dcor(x[!ind1, ] %*% beta_hat[[1]], y[!ind1]) 
  if ((dcor11 + dcor22) >= (dcor12 + dcor21)) {
    return(c(dcor11, dcor22))
  } else {
    return(c(dcor12, dcor21))
  }
  
}

error_rate <- function(est_comp, true_comp) {
  # est_comp: the estimated components, only two clusters
  # true_comp: the true components
  n <- length(est_comp)
  er <- sum(est_comp != true_comp) / n
  if (er > 0.5) {
    #cat("label switch")
    er <- 1 - er
  }
  er
}

error_rate_K <- function(est_comp, true_comp, K) {
  # est_comp: estimated label
  # true_comp: the true label
  # K: the number of clusters
  n <- length(est_comp)
  per_K <- combinat::permn(K) # permutations of c(1:K)
  
  K_prev <- per_K[[1]]
  comp_prev <- K_prev[est_comp]
  for (i in 2:length(per_K)) {
    K_curr <- per_K[[i]]
    comp_curr <- K_curr[est_comp]
    if (sum(true_comp != comp_curr) < sum(true_comp != comp_prev)) {
      K_prev <- K_curr
      comp_prev <- comp_curr
    }
  }
  
  error_rate <- sum(comp_prev != true_comp) / n
  
  list(error_rate = error_rate, pred_K = K_prev, pred_comp = comp_prev)
  
}

error_rate_each_cls <- function(est_comp, true_comp) {
  # For each class, compute the error rate
  # Input:
  # est_comp: the estimated components, only two clusters
  # true_comp: the true components
  n <- length(est_comp)
  n_components <- length(unique(true_comp))
  # compute the number of sample in each cluster
  n_cluster <- numeric(n_components)
  sample_ind_cluster <- vector("list", n_components)
  for (k in seq_len(n_components)) {
    n_cluster[k] <- sum(true_comp == k)
    sample_ind_cluster[[k]] <- which(true_comp == k)
  }
  
  per_K <- combinat::permn(n_components) # permutations of c(1:n_components)
  error_cluster <- rep(Inf, n_components)
  for (i in seq_along(per_K)) {
    K_curr <- per_K[[i]]
    comp_curr <- K_curr[est_comp]
    error_cluster_curr <- rep(NA, n_components)
    for (k in seq_len(n_components)) {
      error_cluster_curr[k] <- sum(comp_curr[sample_ind_cluster[[k]]] != k) / n_cluster[k]
    }
    if (sum(error_cluster_curr) < sum(error_cluster)) {
      error_cluster <- error_cluster_curr
      matched_perm <- K_curr # the matched permutation
    }
  }
  list(error_cluster = error_cluster, matched_perm = matched_perm)
}

switch_est <- function(B_l, comp) {
  # switch the estimated subspace and component in the case of label switching
  # B_l: list with 2 elements
  # comp: vector of predicted component
  B_l <- list(B_l[[2]], B_l[[1]])
  tmp <- comp == 1
  comp[tmp] <- 2
  comp[!tmp] <- 1
  list(B_l = B_l, comp = comp)
}


add_zeros <- function(x, p0) {
  # add p0 zeros to the end of matrix x
  # x: a matrix
  # p0: an integer
  d <- ncol(x)
  rbind(x, matrix(0, nrow = p0, ncol = d))
}


compute_post_ora <- function(X, Y, weight_pi, Gams, Del, etas, FUN, 
                             mu = list(rep(0, ncol(X)), rep(0, ncol(X)))) {
  # compute oracle posterior probability, use true pi, Gam, Delta
  # Input:
  # ===
  # X, Y: data, are n*p matrix and vector with length n
  # weight_pi: a vector consists of pi_1, pi_2
  # Gams: a list, Gams[[1]] = Gam1, Gams[[2]] = Gam2, p*d matrix
  # Del: covariance matrix of the error term
  # etas: d*q matrix
  # FUN: function applied to Y
  # Output:
  # ===
  # a n*2 matrix consists of posterior probabilities
  
  n <- nrow(X)
  d1 <- ncol(Gams[[1]])
  d2 <- ncol(Gams[[2]])
  #mu <- colMeans(X)
  #X_c <- X - matrix(mu, n, ncol(X), byrow = TRUE) # centered predictor
  Fmat <- compute_Fmat(Y, FUN = FUN)
  Fmat_mean <- colMeans(Fmat)
  Fmat_c <- Fmat - matrix(Fmat_mean, nrow(Fmat), ncol(Fmat), byrow = TRUE) # centered function f
  # if (q == 1)
  #   Fmat <- t(Fmat)
  post_prob <- matrix(nrow = n, ncol = 2)
  Del_inv <- solve(Del)
  for (i in seq_len(n)) {
    num_exp <- rep(0, 2) # record the number on the exp
    for (j in seq_len(2)) {
      temp <- matrix(X[i,, drop = FALSE] - mu[[j]] - Fmat_c[i, , drop = FALSE] %*% t(etas[[j]]) %*% t(Gams[[j]]), nrow = 1)
      # to handle the case with too small exp
      num_exp[j] <- 0.5 * temp %*% Del_inv %*% t(temp)
      #post_prob[i, j] <- weight_pi[j] * exp(num_exp[j])
      #cat(i, "-", j, ":", post_prob[i, j], "\n")
    }
    post_prob[i, 1] <- weight_pi[1] / (weight_pi[1] + weight_pi[2]*exp(num_exp[1]-num_exp[2]))
    post_prob[i, 2] <- weight_pi[2] / (weight_pi[1]*exp(num_exp[2]-num_exp[1]) + weight_pi[2])
    #post_prob[i, ] <- post_prob[i, ] / sum(post_prob[i, ])
  }
  post_prob
}

# generalize the function `compute_post_ora` to multiply clusters
compute_post_ora_K <- function(X, Y, weight_pi, Gams, Del, etas, FUN, 
                               mu = NULL) {
  # compute oracle posterior probability, use true pi, Gam, Delta
  # Input:
  # ===
  # X, Y: data, are n*p matrix and vector with length n
  # weight_pi: a vector consists of pi_1, pi_2, ..., pi_K
  # Gams: a list of length K, p*d matrices
  # Del: covariance matrix of the error term
  # etas: a list of d*q matrix
  # FUN: function applied to Y
  # Output:
  # ===
  # a n*K matrix consists of posterior probabilities
  
  n <- nrow(X)
  p <- ncol(X)
  K <- length(weight_pi) # number of clusters
  d_vec <- sapply(Gams, ncol)
  if (is.null(mu)) { 
    mu <- vector("list", K)
    for (k in seq_len(K)) {
      mu[[k]] <- rep(0, p)
    }
  }

  #X_c <- X - matrix(mu, n, ncol(X), byrow = TRUE) # centered predictor
  Fmat <- compute_Fmat(Y, FUN = FUN)
  Fmat_mean <- colMeans(Fmat)
  Fmat_c <- Fmat - matrix(Fmat_mean, nrow(Fmat), ncol(Fmat), byrow = TRUE) # centered function f
  # if (q == 1)
  #   Fmat <- t(Fmat)
  post_prob <- matrix(nrow = n, ncol = K)
  Del_inv <- solve(Del)
  
  for (i in seq_len(n)) {
    num_deno <- rep(0, K) # record the number in the denominator
    for (j in seq_len(K)) {
      for (k in seq_len(K)) {
        if (k == j) {
          temp <- 0
        } else {
          temp <- matrix(X[i,,drop=FALSE] - 0.5*Fmat_c[i,,drop=FALSE]%*%t(Gams[[k]]%*%etas[[k]] + Gams[[j]]%*%etas[[j]]), nrow = 1) # 1*qK
          temp <- temp %*% Del_inv %*% (Gams[[k]]%*%etas[[k]] - Gams[[j]]%*%etas[[j]]) %*% matrix(Fmat_c[i, ], ncol = 1)
        }
        num_deno[j] <- num_deno[j] + weight_pi[j] * exp(temp)
      }
      post_prob[i, j] <- weight_pi[j] / num_deno[j]
      
    }
    
  }
  
  # if NA in gams
  if (sum(is.na(post_prob)) != 0) {
    cat("compute_post_ora_K, update gams, number of NA: ", sum(is.na(post_prob)), "\n")
    post_prob[is.na(post_prob)] <- 0
  }
  
  
  post_prob
}


compute_angle <- function(vector1, vector2) {
  # Compute the dot product of the two vectors
  dot_product <- sum(vector1 * vector2)
  
  # Compute the magnitudes (norms) of the vectors
  magnitude1 <- sqrt(sum(vector1^2))
  magnitude2 <- sqrt(sum(vector2^2))
  
  # Use the dot product and magnitudes to calculate the cosine of the angle
  cosine_theta <- dot_product / (magnitude1 * magnitude2)
  
  # Calculate the angle in radians using the arccos function
  angle_radians <- acos(cosine_theta)
  
  # Convert the angle to degrees
  angle_degrees <- angle_radians * 180 / pi
  
  # if (angle_degrees > 90)
  #   angle_degrees <- 180 - angle_degrees
  return(angle_degrees)
}

# # Example usage:
# vector1 <- c(3, 4, 0)  # Replace with your first vector
# vector2 <- c(1, 2, 0)  # Replace with your second vector
# 
# angle_degrees <- compute_angle(vector1, vector2)
# cat("Angle between the two vectors:", angle_degrees, "degrees\n")
# 
# 




# compute Fmat
compute_Fmat <- function(Y, FUN = NULL) {
  # Y: n-dimensional observation vector for response
  # FUN: the user-specified function f in SEAS-PFC. The default is f(y) = (y, y^2, y^3)
  if (is.null(FUN)) {
    Fmat <- cbind(Y, Y^2, Y^3) # the default function
  } else {
    if (length(FUN(Y[1])) == 1) { # is FUN is 1-d function
      Fmat <- matrix(sapply(Y, FUN), ncol = 1)
    } else {
      Fmat <- t(sapply(Y, FUN)) 
    }
  }
  
  # n <- nrow(Fmat)
  # q <- ncol(Fmat)
  # Fmat_mean <- colMeans(Fmat)
  # Fmat_c <- Fmat - matrix(Fmat_mean, n, q, byrow = TRUE) # centered function matrix
  # lb <- apply(Fmat_c, 2, quantile, 0.1)
  # ub <- apply(Fmat_c, 2, quantile, 0.9)
  # for (i in seq_len(q)) {
  #   Fmat_c[, i] <- sapply(Fmat_c[, i], cut_func, lb[i], ub[i]) # cut extreme value
  # }
  #if (nrow(Fmat == 1)) Fmat <- t(Fmat)
  Fmat
}

# trace function
tr <- function(X) {
  # X: a square matrix
  if (nrow(X) != ncol(X))
    stop("tr function error: X must be a square matrix")
  sum(diag(X))
}


# compute density of multivariate normal distribution
compute_dmvnorm <- function(x, mu, sigma_inv, is_diag = FALSE) {
  # x: a vector
  # mu: mean vector
  # sigma_inv: inverse of covariance matrix
  # is_diag: boolean, true if sigma_inv is a diagnal matrix
  p <- length(mu)
  if (is_diag) {
    det_sigma_inv <- prod(diag(sigma_inv))
  } else {
    det_sigma_inv <- det(sigma_inv)
  }
  
  tmp <- matrix(x - mu, ncol = 1)
  exponent <- -0.5 * t(tmp) %*% sigma_inv %*% tmp
  
  density <- (1 / ((2 * pi)^(p/2)) * sqrt(det_sigma_inv)) * exp(exponent)
  
  return(density)
}

# compute the log of density of multivariate normal distribution
compute_logdmvnorm <- function(x, mu, sigma, perturb = NULL, is_diag = FALSE) {
  if (!is.null(perturb)) {
    diag(sigma) <- diag(sigma) + perturb
  }
  p <- length(mu)
  if (is_diag) {
    det_sigma <- prod(diag(sigma))
  } else {
    det_sigma <- det(sigma)
  }
  
  sigma_inv <- solve(sigma)
  
  tmp <- matrix(x - mu, ncol = 1)
  exponent <- -0.5 * t(tmp) %*% sigma_inv %*% tmp
  
  #log_density <- (1 / ((2 * pi)^(p/2)) * sqrt(det_sigma_inv)) * exp(exponent)
  log_density <- -p/2 * log(2*pi) - 0.5 * log(det_sigma) + exponent
  
  log_density
}

center_Fmat <- function(Fmat, q_lb = 0.1, q_ub = 0.9) {
  Fmat_mean <- colMeans(Fmat)
  Fmat_c <- Fmat - matrix(Fmat_mean, nrow(Fmat), ncol(Fmat), byrow = TRUE) # centered function f
  lb <- apply(Fmat_c, 2, quantile, q_lb)
  ub <- apply(Fmat_c, 2, quantile, q_ub)
  for(i in 1:NCOL(Fmat_c)){
    Fmat_c[, i] <- sapply(Fmat_c[,i], cut_func, lb[i], ub[i])
  }
  Fmat_c
}




compute_mu_K <- function(X, gams) {
  # Input:
  # X: n*p data matrix
  # gams: n*n_components matrix
  # Output:
  # mu: a list with length n_component, each element is a mean vector
  n <- nrow(X) # sample size
  p <- ncol(X) # number of predictors
  n_components <- ncol(gams)
  mu <- vector("list", n_components)
  for (j in seq_len(n_components)) {
    temp <- rep(0, p)
    for (i in seq_len(n)) {
      temp <- temp + gams[i, j] * X[i, ]
    }
    mu[[j]] <- temp / sum(gams[, j])
  }
  mu
}


log_lk <- function(X, Y, B_mats, weight_pi, gams, d_vec = NULL,
                   FUN = NULL, center_F_marg = FALSE) {
  # for estimated beta, compute the loglikelihood
  
  # after the em algorithm finished, we truncate matrices B_k to get beta_k
  # and then estimate the log-likelihood. We estimate it since the p*p covariance
  # matrix Delta is not estimated in the EM algorithm
  # Inputs:
  # ======
  # X: n*p design matrix
  # Y: length n vector
  # B_mats: a list of p*q matrix
  # gams: n*n_components
  # d_vec: true dimensions of Gamma1, Gamma2, ..., GammaK
  # Outputs:
  # =======
  # the log_lk
  if (is.null(d_vec)) {
    stop("log_lk: true dimension d_vec is missing")
  }
  
  Fmat <- compute_Fmat(Y, FUN) 
  q <- ncol(Fmat)
  if (q == nrow(X))
    Fmat <- t(Fmat)
  q <- ncol(Fmat)
  n_components <- length(d_vec) # number of clusters
  
  beta_hat_mats <- vector("list", length = n_components)
  for (k in seq_len(n_components)) {
    d_k <- d_vec[k]
    beta_hat_mats[[k]] <- svd(B_mats[[k]])$u[, 1:d_k, drop = FALSE]
  }
  
  # now we compute the loglikelihood of Z, we just do one more e-step using
  # function `e_step_K`
  params <- list(B_mats = beta_hat_mats, gams = gams, weight_pi = weight_pi)
  loglk <- e_step_K(X, Fmat, params, center_F_marg = center_F_marg, compute_loglk = TRUE)$loglk
  
  loglk
}


log_lk2 <- function(X, Y, B_mats, weight_pi, gams, d_vec = NULL,
                   FUN = NULL, center_F_marg = FALSE) {
  # for estimated beta, compute the loglikelihood
  
  # after the em algorithm finished, we truncate matrices B_k to get beta_k
  # and then estimate the log-likelihood. We estimate it since the p*p covariance
  # matrix Delta is not estimated in the EM algorithm
  # Inputs:
  # ======
  # X: n*p design matrix
  # Y: length n vector
  # B_mats: a list of p*q matrix
  # gams: n*n_components
  # d_vec: true dimensions of Gamma1, Gamma2, ..., GammaK
  # Outputs:
  # =======
  # the log_lk
  if (is.null(d_vec)) {
    stop("log_lk: true dimension d_vec is missing")
  }
  
  Fmat <- compute_Fmat(Y, FUN) 
  q <- ncol(Fmat)
  if (q == nrow(X))
    Fmat <- t(Fmat)
  q <- ncol(Fmat)
  n_components <- length(d_vec) # number of clusters
  
  beta_hat_mats <- vector("list", length = n_components)
  for (k in seq_len(n_components)) {
    d_k <- d_vec[k]
    beta_hat_mats[[k]] <- svd(B_mats[[k]])$u[, 1:d_k, drop = FALSE]
  }
  beta_hat <- do.call(cbind, beta_hat_mats) # p * sum(d_vec)
  # compute the orthogonal complement
  beta_0_hat <- orthogonal_complement(beta_hat) 
  beta_0_hat_mats <- vector("list", length = n_components)
  for (k in seq_len(n_components)) {
    beta_0_hat_mats[[k]] <- beta_0_hat
  }
  
  # now we compute the loglikelihood of Z, we just do one more e-step using
  # function `e_step_K`
  params1 <- list(B_mats = beta_hat_mats, gams = gams, weight_pi = weight_pi)
  loglk1 <- e_step_K(X, Fmat, params1, center_F_marg = center_F_marg, compute_loglk = TRUE)$loglk
  
  params2 <- list(B_mats = beta_0_hat_mats, gams = gams, weight_pi = weight_pi)
  loglk2 <- e_step_K(X, Fmat, params2, center_F_marg = center_F_marg, compute_loglk = TRUE)$loglk
  
  loglk1 + loglk2
}


orthogonal_complement <- function(A) {
  # get the orthogonal subspace of column space of A
  n <- nrow(A)
  # Perform SVD decomposition
  svd_A <- svd(A, nu = n)
  
  # The columns of U corresponding to zero singular values (in D) span the orthogonal complement
  rank_A <- sum(svd_A$d > 1e-6)  # Tolerance to detect non-zero singular values
  if (rank_A == n) {
    complement <- matrix(0, nrow = n, ncol = 1)
  } else {
    complement <- svd_A$u[, (rank_A + 1):n, drop = FALSE]  # Drop=FALSE keeps matrix structure
  }

  complement
}



est_gam <- function(X, Y, params, FUN = NULL) {
  # compute the estimated eta
  # for new input Xnew and Ynew, we can predict their cluster
  # this function is similar to the e_step
  # input:
  # X, Y: the variable matrix and response 
  # params: a list of parameters used to compute eta
  #        B_mats: a list of p*q matrix, length equals to number of clusters
  #        A_mats: a list of qK * q matrix, length equals to number of clusters
  #        Delta_inv: qK * qK matrix
  #        weight_pi: length n_components vector
  #        mu: a list of cluster means
  # FUN: a function to compute Fmat
  if (is.null(FUN)) {
    stop("est_gam: FUN is NULL!")
  }
  # unpack params
  B_mats <- params$B_mats
  weight_pi <- params$weight_pi
  A_mats <- params$A_mats
  Delta_inv <- params$Delta_inv
  mu <- params$mu
  K <- length(weight_pi) # number of clusters
  n <- nrow(X) # number of samples
  p <- ncol(X) # number of variables
  
  if (n != length(Y)) {
    stop("X and Y do not have the same number of samples!")
  }
  
  # check if FUN is valid
  Fmat <- compute_Fmat(Y, FUN) 
  q <- ncol(Fmat)
  if (q != ncol(B_mats[[1]])) {
    stop("est_gam: wrong FUN function! The FUN should have the same dimenion of training FUN!")
  }
  
  gams <- matrix(nrow = n, ncol = K) # initialize gam matrix
  # compute centered X 
  X_c <- vector("list", K) 
  for (j in seq_len(K)) {
    X_c[[j]] <- X - matrix(mu[[j]], n, p, byrow = TRUE)
  }
  # compute centered Fmat
  Fmat_c <- center_Fmat(Fmat) # Fmat is centered marginally 
  
  # compute Z mats
  Z_mats <- vector("list", length = K)
  for (j in seq_len(K)) {
    Z_mats[[j]] <- X_c[[j]] %*% B_mats[[j]] # n * q
  }
  Z <- do.call(cbind, Z_mats) # n * qK
  # compute gams
  for (i in seq_len(n)) {
    num_deno <- rep(0, K) # record the number in the denominator
    for (j in seq_len(K)) {
      for (k in seq_len(K)) {
        if (k == j) {
          temp <- 0
        } else {
          temp <- matrix(Z[i,,drop=FALSE] - 0.5*matrix(Fmat_c[i, ], nrow = 1)%*%t(A_mats[[k]] + A_mats[[j]]), nrow = 1) # 1*qK
          temp <- temp %*% Delta_inv %*% (A_mats[[k]] - A_mats[[j]]) %*% matrix(Fmat_c[i, ], ncol = 1)
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
  
  gams
}













