#' Produce sufficient statistics (lambda0) given the common coefficients (rho0)
#'
#' @param rho lagged dependent variable coefficients
#' @param alpha external variables coefficients
#' @param N cross-sectional dimension of the data
#' @param T time dimension of the data
#' @param n_alpha number of external variables
#' @param Y_mat dependent variable matrix
#' @param X_mat lagged dependent variable matrix
#' @param W cross-sectionally invariant variables - not used now
#' @param Z_mat external variable matrix
#'
#' @importFrom pracma mldivide
get_lambda0 <- function(rho, alpha = rep(0, n_alpha), N, T, n_alpha, Y_mat,
                        X_mat, W, Z_mat) {
  Y_tilde <- matrix(NA, nrow = T, ncol = N)
  for (t in 1:T) {
    if (!is.null(Z_mat)) {
      Z_fitted_t <- rep(0, N)
      for (a in 1:n_alpha) {
        Z_fitted_t <- Z_fitted_t + c(as.matrix(alpha[a]) %*%
          t(as.matrix(Z_mat[[a]][t, ])))
      }
      Y_tilde[t, ] <- Y_mat[t, ] - rho %*% X_mat[t, ] - Z_fitted_t
    } else {
      Y_tilde[t, ] <- Y_mat[t, ] - rho %*% X_mat[t, ]
    }
  }
  lambda_0 <- mldivide((t(W) %*% W), (t(W) %*% Y_tilde))
  return(lambda_0)
}

#' Produce variance of the shocks estimated using GMM residues (sigma2_0)
#' given the common coefficients (rho0)
#'
#' @param rho lagged dependent variable coefficients
#' @param alpha external variables coefficients
#' @param common_par_method method for estimating common parameters
#' @param X_star auxiliary matrix for OFD transformation
#' @param Y_star auxiliary matrix for OFD transformation
#' @param Z_star auxiliary matrix for OFD transformation
#' @param X_mat lagged dependent variable matrix
#' @param Y_mat dependent variable matrix
#' @param Z_mat external variable matrix
#' @param n_alpha number of external variables
get_sigma2 <- function(rho, alpha = 0, common_par_method, X_star, Y_star, Z_star,
                       X_mat, Y_mat, Z_mat, n_alpha) {
  if (common_par_method == "GMM.ABover") {
    if (!is.null(Z_star)) {
      Z_fitted <- rep(0, nrow(Y_star) * ncol(Y_star))
      for (a in 1:n_alpha) {
        Z_fitted <- Z_fitted + as.numeric(alpha[a]) %*% c(as.matrix(Z_star[[a]]))
      }
      sigma2_0 <- var(t(c(Y_star) - rho %*% c(X_star) - Z_fitted))
    } else {
      sigma2_0 <- var(t(c(Y_star) - rho %*% c(X_star)))
    }
  } else {
    if (!is.null(Z_star)) {
      Z_fitted <- rep(0, nrow(Y_mat) * ncol(Y_mat))
      for (a in 1:n_alpha) {
        Z_fitted <- Z_fitted + as.numeric(alpha[a]) %*% c(as.matrix(Z_mat[[a]]))
      }
      sigma2_0 <- var(t(c(Y_mat) - rho %*% c(X_mat) - Z_fitted))
    } else {
      sigma2_0 <- var(t(c(Y_mat) - rho %*% c(X_mat)))
    }
  }
  return(sigma2_0)
}

#' Produce posterior means of lambda's for the parametric GMM implementation
#' given autoregressive coefficient (rho)
#'
#' @param rho lagged dependent variable coefficients
#' @param alpha external variables coefficients
#' @param optim_method optimization method
#' @param init initial values for the optimization routine
#' @param n_lambda number of columns in W; currently always set to 1
#' @param n_alpha number of external variables
#' @param X_mat lagged dependent variable matrix
#' @param Y_mat dependent variable matrix
#' @param Z_mat external variable matrix
#' @param W cross-sectionally invariant variables - not used now
#' @param T time dimension of the data
#' @param N cross-sectional dimension of the data
#' @param aux_Y0 auxiliary matrix with initial observations of the dependent variable 
#' @param common_par_method method for estimating common parameters
#' @param X_star auxiliary matrix for OFD transformation
#' @param Y_star auxiliary matrix for OFD transformation
#' @param Z_star auxiliary matrix for OFD transformation
#'
#' @importFrom minqa bobyqa
GMM_parametric <- function(rho, alpha = 0, optim_method, init, n_lambda, n_alpha,
                           X_mat, Y_mat, Z_mat, W, T, N, aux_Y0, common_par_method,
                           X_star, Y_star, Z_star) {
  rho_GMMpar <- c(rho)
  alpha_GMMpar <- c(alpha)
  lambda_GMMpar <- get_lambda0(
    rho_GMMpar, alpha_GMMpar, N, T, n_alpha, Y_mat,
    X_mat, W, Z_mat
  )
  sigma2_GMMpar <- c(get_sigma2(
    rho_GMMpar, alpha_GMMpar, common_par_method, X_star,
    Y_star, Z_star, X_mat, Y_mat, Z_mat, n_alpha
  ))
  if (optim_method == "gradient") {
    minimise_GMM <- try(minimise_GMM <- optim(init[(n_alpha + 3):length(init)],
      loglikelihood_GMM,
      rho_GMMpar = rho_GMMpar,
      alpha_GMMpar = alpha_GMMpar,
      sigma2_GMMpar = sigma2_GMMpar,
      n_alpha = n_alpha, X_mat = X_mat,
      Y_mat = Y_mat, Z_mat = Z_mat,
      W = W, T = T, N = N,
      aux_Y0 = aux_Y0
    ))
    if ("try-error" %in% class(minimise_GMM)) {
      print("Gradient optimisation failed. Quadratic approximation was used.")
      minimise_GMM <- bobyqa(init[(n_alpha + 3):length(init)], loglikelihood_GMM,
        rho_GMMpar = rho_GMMpar,
        alpha_GMMpar = alpha_GMMpar,
        sigma2_GMMpar = sigma2_GMMpar,
        n_alpha = n_alpha, X_mat = X_mat,
        Y_mat = Y_mat, Z_mat = Z_mat,
        W = W, T = T, N = N,
        aux_Y0 = aux_Y0
      )
    }
  } else if (optim_method == "quadratic") {
    minimise_GMM <- bobyqa(init[(n_alpha + 3):length(init)], loglikelihood_GMM,
      rho_GMMpar = rho_GMMpar,
      alpha_GMMpar = alpha_GMMpar,
      sigma2_GMMpar = sigma2_GMMpar,
      n_alpha = n_alpha, X_mat = X_mat,
      Y_mat = Y_mat, Z_mat = Z_mat,
      W = W, T = T, N = N,
      aux_Y0 = aux_Y0
    )
  } else if (optim_method == "annealing") {
    minimise_GMM <- optim(init[(n_alpha + 3):length(init)], loglikelihood_GMM,
      rho_GMMpar = rho_GMMpar,
      alpha_GMMpar = alpha_GMMpar,
      sigma2_GMMpar = sigma2_GMMpar,
      n_alpha = n_alpha, X_mat = X_mat,
      Y_mat = Y_mat, Z_mat = Z_mat,
      W = W, T = T, N = N,
      aux_Y0 = aux_Y0,
      method = "SANN"
    )
  }
  param_GMM <- minimise_GMM$par
  mmu_GMMpar <- matrix(param_GMM[1:(length(param_GMM) - n_lambda)], nrow = 2)
  ww2_lambda_GMMpar <- param_GMM[(length(param_GMM) - n_lambda + 1):
  length(param_GMM)]
  mean_lambda_GMMpar <- post_mean_lambda_par(
    lambda_GMMpar, sigma2_GMMpar,
    mmu_GMMpar, ww2_lambda_GMMpar, W,
    aux_Y0
  )
  return(mean_lambda_GMMpar)
}

#' Obtain 2D kernel density estimates given sufficient statistics for lambdas
#' and the initial data Y0
#'
#' @param lambdas sufficient statistics for the intercept term
#' @param sigma2 variance of the shocks
#' @param dens_grid grid over which the density is to be computed
#' @param N cross-sectional dimension of the data
#' @param T time dimension of the data
#' @param Y0 initial observations of the dependent variable 
#'
#' @importFrom pracma interp2 ones
get_kernel <- function(lambdas, sigma2, dens_grid, N, T, Y0) {
  zY0 <- cbind(t(lambdas), as.matrix(Y0))
  max <- c(max(zY0[, 1]), max(zY0[, 2]))
  min <- c(min(zY0[, 1]), min(zY0[, 2]))
  range <- max - min
  max_xy <- max + range / 100
  min_xy <- min - range / 100
  kernel <- kde2D(zY0, n = dens_grid, limits = c(
    min_xy[1], max_xy[1],
    min_xy[2], max_xy[2]
  ))
  q_zy <- kernel$density
  xx <- kernel$X
  yy <- kernel$Y
  kernelY0 <- kde(Y0, n = dens_grid, MIN = min(yy[, 1]), MAX = max(yy[, 1]))
  q_y <- kernelY0[2, ]
  n_qz <- ncol(q_zy)
  qq <- q_zy / (q_y %*% ones(1, n_qz))
  i_aux <- rowSums(ones(N, 1) %*% xx[1, ] < t(lambdas) %*% ones(1, n_qz))
  q_aux <- interp2(xx[1, ], t(yy[, 1]), qq, cbind(
    as.matrix(xx[1, i_aux]),
    t(lambdas), as.matrix(xx[1, i_aux + 1])
  ),
  as.matrix(Y0) %*% ones(1, 3),
  method = "linear"
  )
  q_aux <- matrix(q_aux, ncol = 3)
  postmean_lambda <- lambdas + t((((q_aux[, 3] - q_aux[, 1]) /
    (xx[1, i_aux + 1] - xx[1, i_aux])) / (q_aux[, 2])) %*%
    as.matrix(sigma2) / T)
  return(postmean_lambda)
}

#' Produce (negative) log marginal likelihood for QMLE with correlated random
#' coefficients
#'
#' @param param vectores of parameters to optimize over
#' @param n_alpha number of external variables
#' @param X_mat lagged dependent variable matrix
#' @param Y_mat dependent variable matrix
#' @param Z_mat external variable matrix
#' @param W cross-sectionally invariant variables - not used now
#' @param T time dimension of the data
#' @param N cross-sectional dimension of the data
#' @param aux_Y0 auxiliary matrix with initial observations of the dependent variable 
#'
#' @importFrom pracma mldivide mrdivide Diag
loglikelihood_QMLE <- function(param, n_alpha, X_mat, Y_mat, Z_mat, W, T, N,
                               aux_Y0) {
  # construct parameters vector
  n_param <- length(param)
  n_lambda <- (n_param - 2 - n_alpha) / 3
  rho_0 <- param[1]
  if (n_alpha == 0) {
    alpha_0 <- 0
  } else {
    alpha_0 <- param[2:(n_alpha + 1)]
  }
  sigma2_0 <- param[n_alpha + 2]
  mmu <- matrix(param[(n_alpha + 3):(length(param) - n_lambda)], nrow = 2)
  ww2_lambda <- param[(length(param) - n_lambda + 1):length(param)]
  # standardise data
  X_stand <- X_mat / sqrt(sigma2_0)
  Y_stand <- Y_mat / sqrt(sigma2_0)
  if (!is.null(Z_mat)) {
    Z_stand <- vector("list", n_alpha)
    for (a in 1:n_alpha) {
      Z_stand[[a]] <- Z_mat[[a]] / sqrt(sigma2_0)
    }
  }
  W_stand <- W / sqrt(sigma2_0)
  # construct lambda_0
  Y_tilde_stand <- matrix(NA, nrow = T, ncol = N)
  for (t in 1:T) {
    if (!is.null(Z_mat)) {
      Z_fit_t_stand <- rep(0, N)
      for (a in 1:n_alpha) {
        Z_fit_t_stand <- Z_fit_t_stand + c(as.matrix(alpha_0[a]) %*%
          t(as.matrix(Z_stand[[a]][t, ])))
      }
      Y_tilde_stand[t, ] <- Y_stand[t, ] - rho_0 %*% X_stand[t, ] - Z_fit_t_stand
    } else {
      Y_tilde_stand[t, ] <- Y_stand[t, ] - rho_0 %*% X_stand[t, ]
    }
  }
  lambda_0 <- mldivide((t(W_stand) %*% W_stand), (t(W_stand) %*% Y_tilde_stand))
  # construct sum_post
  lambda_n <- t(W_stand) %*% W_stand %*% lambda_0 +
    Diag(1 / ww2_lambda) %*% t(mmu) %*% aux_Y0
  lambda_d <- t(W_stand) %*% W_stand + Diag(1 / ww2_lambda)
  sum_post <- 0
  for (i in 1:N) {
    sum_post <- sum_post + mrdivide(t(lambda_n[, i]), lambda_d) %*% lambda_n[, i]
  }
  # calculate loglikelihood
  L <- suppressWarnings(
    N * T * log(sigma2_0) + N * sum(log(ww2_lambda)) + N * log(det(lambda_d)) +
      sum(c(Y_tilde_stand) ^ 2) + sum(colSums((t(mmu) %*% aux_Y0) ^ 2) / ww2_lambda) -
      sum_post
  )
  return(L)
}

#' Produce negative log-likelihood in the GMM case
#'
#' @param theta vector of homogeneous parameters
#' @param rho_GMMpar lagged dependent variables coefficient estimates from the GMM
#' @param alpha_GMMpar external variables coefficient estimates from the GMM
#' @param sigma2_GMMpar variance of the shocks estimated using GMM residuals
#' @param n_alpha number of external variables
#' @param X_mat lagged dependent variable matrix
#' @param Y_mat dependent variable matrix
#' @param Z_mat external variable matrix
#' @param W cross-sectionally invariant variables - not used now
#' @param T time dimension of the data
#' @param N cross-sectional dimension of the data
#' @param aux_Y0 auxiliary matrix with initial observations of the dependent variable 
loglikelihood_GMM <- function(theta, rho_GMMpar, alpha_GMMpar, sigma2_GMMpar,
                              n_alpha, X_mat, Y_mat, Z_mat, W, T, N, aux_Y0) {
  if (n_alpha > 0) {
    param <- c(rho_GMMpar, alpha_GMMpar, sigma2_GMMpar, theta)
  }
  else {
    param <- c(rho_GMMpar, sigma2_GMMpar, theta)
  }
  L <- loglikelihood_QMLE(param, n_alpha, X_mat, Y_mat, Z_mat, W, T, N, aux_Y0)
  return(L)
}

#' Provide posterior means of lambda_i's based on the Parametric Posterior Mean
#' estimator with correlated random coefficients
#'
#' @param lambda0 initial estimate of lambdas
#' @param sigma2 variance of the shocks
#' @param mmu auxiliary result (mean)
#' @param ww2_lambda auxiliary result (lambda times ww2)
#' @param W cross-sectionally invariant variables - not used now
#' @param aux_Y0 auxiliary matrix with initial observations of the dependent variable 
#'
#' @importFrom pracma mldivide Diag
post_mean_lambda_par <- function(lambda0, sigma2, mmu, ww2_lambda, W, aux_Y0) {
  n_lambda <- ncol(mmu)
  W <- W / sigma2
  m_lambda <- mldivide(
    (t(W) %*% W + Diag(1 / ww2_lambda)),
    t(W) %*% W %*% lambda0 + Diag(1 / ww2_lambda) %*% t(mmu) %*% aux_Y0
  )
  return(m_lambda)
}
