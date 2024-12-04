# ------------------ revised functions from 'msda' package ---------------------- #
# Revise 'msda' function to accommodate other forms of M and U matrices.
# Some inputs are similar to the ones in 'seas' function. Please refer to 'msda' package documentation for more details of the arguments in 'msda' function.
# 
# Outputs:
# ========
# lambda: The tuning parameter sequence.
# theta: The list of estimated matrix.
# M: The M matrix from samples. It is the sample covariance of predictor for SEAS-SIR, SEAS-Intra, and SEAS-PFC.
# U: The U matrix from samples, which depends on the argument type.
# rank: The list of estimated rank for each matrix.
########################## following paras are for mixpfc ######################
# gam: n dimensional vector, which is the j-th column of gams (n*n_components matrix)
msda <- function(x, y, yclass=NULL, categorical=FALSE, H=5, type='sir', FUN = NULL, 
                 gam = NULL, center_F_marg = center_F_marg,
                 lambda.factor=NULL, nlambda=100, lambda=NULL, dfmax=NULL, pmax=NULL,
                 pf=NULL, M = NULL, U = NULL, nobs=NULL, nclass=NULL, eps=1e-04,
                 maxit=1e+06, sml=1e-06, verbose=FALSE, perturb=NULL){
  if(is.null(M) || is.null(U)){ # Generate M and U matrices
    if(missing(x) || missing(y)) stop("Missing x or y.")
    if(is.data.frame(x)) x <- as.matrix(x)
    if(is.null(yclass)){
      if(categorical == FALSE){
        ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
        yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
      }
      else if(categorical == TRUE){
        yclass <- y
      }
    }
    if(any(table(yclass) < 5)) warning(sprintf("The sample size of class %d is less than 5\n", which(table(yclass) < 5)))
    nclass <- as.integer(length(unique(yclass)))
    MU_out <- MU(x, y, yclass, type, FUN, gam = gam, center_F_marg = center_F_marg)
    M <- MU_out$M
    U <- MU_out$U
    nobs <- as.integer(dim(x)[1])
    nvars <- as.integer(dim(x)[2])
  }
  else{
    if(is.null(nobs)) stop("Missing nobs.")
    if(is.null(nclass)) stop("Missing nclass.")
    nvars <- NCOL(M)
  }
  
  if(is.null(lambda.factor)) lambda.factor <- ifelse((nobs - nclass)<=nvars, 0.2, 1e-03)
  if(is.null(dfmax)) dfmax <- nobs
  if(is.null(pmax)) pmax <- min(dfmax*2 + 20, nvars)
  if(is.null(pf)) pf <- rep(1, nvars)
  if (!is.null(perturb)) 
    diag(M) <- diag(M) + perturb
  H <- as.integer(dim(U)[2])
  ## parameter setup
  if (length(pf) != nvars) 
    stop("The size of penalty factor must be same as the number of input variables")
  maxit <- as.integer(maxit)
  verbose <- as.integer(verbose)
  sml <- as.double(sml)
  pf <- as.double(pf)
  eps <- as.double(eps)
  dfmax <- as.integer(dfmax)
  pmax <- as.integer(pmax)
  ## lambda setup
  nlam <- as.integer(nlambda)
  if (is.null(lambda)) {
    if (lambda.factor >= 1)
      stop("lambda.factor should be less than 1")
    flmin <- as.double(lambda.factor)
    ulam <- double(1)  #ulam=0 if lambda is missing
  } else {
    # flmin=1 if user define lambda
    flmin <- as.double(1)
    if (any(lambda < 0))
      stop("lambdas should be non-negative")
    ulam <- as.double(rev(sort(lambda)))
    nlam <- as.integer(length(lambda))
  }
  ## call Fortran core
  fit <- .Fortran("msda", obj = double(nlam), H, nvars, as.double(M), as.double(t(U)), pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, sml, verbose, nalam = integer(1), theta = double(pmax * H * nlam), itheta = integer(pmax), ntheta = integer(nlam),alam = double(nlam), npass = integer(1), jerr = integer(1))
  ## output
  outlist <- formatoutput(fit, maxit, pmax, nvars, H)
  rank <- rep(NA_integer_, length(outlist$theta))
  for (i in 1:length(outlist$theta)){
    if(!is.null(outlist$theta[[i]])){
      rank[i] <- rank_func(outlist$theta[[i]], thrd = 1e-3)
    }
  }
  if(is.null(lambda))
    outlist$lambda <- lamfix(outlist$lambda)
  outlist <- list(lambda = outlist$lambda, theta = outlist$theta, M = M, U = U, rank = rank)
  class(outlist) <- c("msda")
  outlist
}

# Revise 'cv.msda' function to accommodate other forms of M and U matrices. We also add the optional argument 'fold' to pass the user-specified folds index.
# Some inputs are similar to the ones in 'cv.seas' function. Please refer to 'msda' package documentation for more details of arguments used in 'cv.msda' function.
# 
# Outputs:
# ========
# beta: The optimal estimated matrix.
# id: The index of the optimal tuning parameter.
# lambda: The lambda sequence.
# lam_max: The optimal tuning parameter.
# rank: The rank of the optimal estimated matrix.
# M: The M matrix based on the full data. It is the sample covariance of predictor for SEAS-SIR, SEAS-Intra, and SEAS-PFC.
# U: The U matrix based on the full data, which depends on the argument type.
# M_fold: The M matrix list based on each cross-validation data fold.
# U_fold: The U matrix list based on each cross-validation data fold.
########################## following paras are for mixpfc ######################
# gam: n dimensional vector, which is the j-th column of gams (n*n_components matrix)
cv.msda <- function(x, y, yclass=NULL, categorical=FALSE, H=5, type='sir', 
                    gam = NULL, center_F_marg = FALSE,
                    lambda.factor=NULL, nlambda=100, nfolds=5, foldid = NULL, 
                    lambda = NULL, FUN = NULL, maxit = 1e3, plot = FALSE){
  if(is.data.frame(x)) x <- as.matrix(x)
  if(is.null(yclass)){
    if(categorical == FALSE){
      ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
      yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
    }
    else if(categorical == TRUE){
      yclass <- y
    }
  }
  if(any(table(yclass) < 5)) warning(sprintf("The sample size of class %d is less than 5\n", which(table(yclass) < 5)))
  nobs <- nrow(x)
  nvars <- ncol(x)
  nclass <- length(unique(yclass))
  if(is.null(lambda.factor)) lambda.factor <- ifelse((nobs - nclass)<=nvars, 0.2, 1e-03)
  # Fit the model on the full data, obtain the lambda sequence.
  fit <- msda(x, y, yclass = yclass, type = type, gam = gam, center_F_marg = center_F_marg,
              lambda.factor = lambda.factor, 
              nlambda = nlambda, lambda = lambda, FUN = FUN, maxit=maxit)
  lambda <- fit$lambda
  beta_l <- fit$theta
  M <- fit$M
  U <- fit$U
  rank_l <- fit$rank
  beta_l <- cut_mat(beta_l, 1e-3, rank_l)
  
  # Cross-validation
  if(is.null(foldid)){
    ord <- order(y)
    y <- y[ord]
    yclass <- yclass[ord]
    x <- x[ord,]
    gam <- gam[ord] # gam 
    count <- as.numeric(table(yclass))
    foldid <- c()
    for(cnt in count){
      foldid <- c(foldid, sample(rep(seq(nfolds), length = cnt)))
    }
  }
  else{
    nfolds <- length(unique(foldid))
  }

  cv_out <- lapply(1:nfolds, function(k){
    x_train <- x[foldid!=k,,drop=FALSE]
    x_val <- x[foldid==k,,drop=FALSE]
    y_train <- y[foldid!=k]
    y_val <- y[foldid==k]
    yclass_train <- yclass[foldid!=k]
    #####
    gam_train <- gam[foldid!=k]
    gam_val <- gam[foldid==k]
    #####
    fit_fold <- msda(x_train, y_train, yclass_train, type = type, 
                     gam = gam_train, center_F_marg = center_F_marg,
                     lambda.factor=lambda.factor, nlambda=nlambda, lambda = lambda, 
                     FUN = FUN, maxit=maxit)
    M_fold <- fit_fold$M
    U_fold <- fit_fold$U
    beta_fold <- fit_fold$theta
    rank_fold <- fit_fold$rank
    beta_fold <- cut_mat(beta_fold, 1e-3, rank_fold)
    
    # return evaluation of each fold
    eval_fold <- eval_dc(beta_fold, x_val, y_val, type = type, gam = gam_val)
    if(length(eval_fold) != length(lambda)){
      eval_fold <- c(eval_fold, rep(NA, length(lambda) - length(eval_fold)))
    }
    list(eval = eval_fold, M = M_fold, U = U_fold)
  })
  
  eval_all <- do.call(rbind, lapply(cv_out, "[[", 1))
  M_fold <- lapply(cv_out, "[[", 2)
  U_fold <- lapply(cv_out, "[[", 3)
  if(is.vector(eval_all)){
   eval_all <- t(as.matrix(eval_all))
  }

  ## No matrix is converged in any fold
  if(all(is.na(eval_all))) return(NULL)
  
  cvm <- colMeans(eval_all, na.rm = TRUE)
  # The optimal lambda1
  id_max <- which.max(cvm)
  lam_max <- lambda[id_max]
  beta <- as.matrix(beta_l[[id_max]])
  
  # Recalculate the rank
  rank <- rank_func(beta, thrd = 1e-3)
  
  if(plot){ # If TRUE, plot the cv evaluation for each tuning parameter.
    dat <- data.frame(x = 1:length(cvm), y = cvm)
    g <- ggplot(dat, aes(x = x, y = y))+
      geom_point(size = 1)+
      xlab("")+
      ylab("Distance correlation")+
      theme_bw()
    print(g)
  }
  
  list(beta = beta, beta_l = beta_l, id = id_max, lambda = lambda, lam_max = lam_max, rank = rank, M = M, U = U, M_fold = M_fold, U_fold = U_fold)
}
