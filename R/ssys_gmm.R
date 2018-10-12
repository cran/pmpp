#' Suboptimal multi-step System-GMM estimator for AR(1) panel data model
#'
#' @author Michal Oleszak
#'
#' @description Computes an enhanced version of the Blundell-Bond (System-GMM)
#'              estimator for panel data by means of replacing the standard
#'              GMM-weighting matrix by its sub-optimal version, thus increasing
#'              estimator's efficiency.
#'
#' @references Youssef, A. and Abonazel, M. (2015). Alternative GMM estimators
#'             for first-order autoregressive panel model: An improving
#'             efficiency approach. MPRA Paper No. 68674; Forthcoming in:
#'             Communications in Statistics - Simulation and Computation,
#'             \url{https://mpra.ub.uni-muenchen.de/68674/1/MPRA_paper_68674.pdf}
#'
#' @param Y     matrix of size (T x N) with the dependent variable
#' @param model one of: onestep, twosteps, threesteps; more steps should
#'              increase efficiency, but might be computationally infeasible
#'              (a singular matrix needs to be inverted); if this is the case,
#'              generalised inverse is used
#'
#' @return The estimated value of the auto-regressive parameter.
#'
#' @importFrom plm pgmm
#' @importFrom Matrix bdiag solve
#' @importFrom pracma ones
#' @importFrom MASS ginv
#' @importFrom utils head tail
#' @export
ssys_gmm <- function(Y, model = c("onestep", "twosteps", "threesteps")) {
  N <- ncol(Y)
  T <- nrow(Y)
  # tranform data to a format accepted by pgmm()
  unit <- rep(1:N, each = T)
  time <- rep(1:T, N)
  yplm <- c(Y)
  plmdata <- cbind(unit, time, yplm)
  # obtain variance ratio 'r'
  DIF <- pgmm(yplm ~ lag(yplm, 1) | lag(yplm, 2:99),
    data = as.data.frame(plmdata), effect = "individual",
    method = "onestep", transformation = "d"
  )
  SYS <- pgmm(yplm ~ lag(yplm, 1) | lag(yplm, 2:99),
    data = as.data.frame(plmdata), effect = "individual",
    method = "onestep", transformation = "ld"
  )
  delta_u_i_hat <- DIF$residuals
  sysres <- SYS$residuals
  delta_u_i_tilde <- lapply(sysres, function(x) head(x, (T - 2)))
  u_i_tilde <- lapply(sysres, function(x) tail(x, (T - 1)))
  sigma2eps <- (sum(unlist(lapply(delta_u_i_hat, function(x) t(x) %*% x)))) /
    (2 * N * (T - 2))
  sigma2gamma <- (sum(unlist(lapply(u_i_tilde, function(x) t(x) %*% x)) -
    (unlist(lapply(delta_u_i_tilde, function(x) t(x) %*% x)) / 2))) /
    (N * (T - 2))
  r <- sigma2gamma / sigma2eps
  # construct moment condition matrices
  L_i <- list()
  D_i <- list()
  D_i_list <- rep(list(vector("list", T - 2)), N)
  S_i <- list()
  Y_diff <- matrix(NA, ncol = N, nrow = T - 2)
  for (t in 2:(T - 1)) {
    Y_diff[t - 1, ] <- Y[t, ] - Y[t - 1, ]
  }
  for (i in 1:N) {
    L_i[[i]] <- diag(Y_diff[, i])
    for (t in 1:(T - 2)) {
      D_i_list[[i]][[t]] <- matrix(c(Y[1:t, i]), nrow = 1)
    }
    D_i[[i]] <- bdiag(D_i_list[[i]])
    S_i[[i]] <- bdiag(D_i[[i]], L_i[[i]])
  }
  F <- matrix(0, nrow = T - 2, ncol = T - 1)
  for (i in 1:nrow(F)) {
    for (j in 1:ncol(F)) {
      if (i == j) {
        F[i, j] <- -1
      }
      if (j == i + 1) {
        F[i, j] <- 1
      }
    }
  }
  A <- diag(T - 2) + r * ones(T - 2)
  G <- bdiag(F %*% t(F), A)
  Z1 <- matrix(0, nrow = ncol(S_i[[1]]), ncol = ncol(S_i[[1]]))
  for (i in 1:N) {
    Z1 <- Z1 + t(as.matrix(S_i[[i]])) %*% G %*% S_i[[i]]
  }
  Z1 <- solve(Z1)
  S <- S_i[[1]]
  for (i in 2:N) {
    S <- rbind(S, S_i[[i]])
  }
  Y_diff_temp <- matrix(NA, ncol = N, nrow = T - 2)
  for (t in 3:(T)) {
    Y_diff_temp[t - 2, ] <- Y[t, ] - Y[t - 1, ]
  }
  y_stacked <- c(c(Y_diff_temp), c(Y[-c(1:2), ]))
  y_stacked_lag <- c(c(Y_diff), c(Y[-c(1, T), ]))
  rho1 <- solve(t(y_stacked_lag) %*% S %*% Z1 %*% t(as.matrix(S)) %*% y_stacked_lag) %*%
    (t(y_stacked_lag) %*% S %*% Z1 %*% t(as.matrix(S)) %*% y_stacked)
  if (model == "onestep") {
    return(as.numeric(rho1))
  }

  # step 2
  fitted1lev <- matrix(NA, nrow = T - 2, ncol = N)
  fitted1diff <- matrix(NA, nrow = T - 2, ncol = N)
  for (t in 1:(T - 2)) {
    fitted1lev[t, ] <- Y[t + 1, ] * as.numeric(rho1)
    fitted1diff[t, ] <- Y_diff_temp[t, ] * as.numeric(rho1)
  }
  res1lev <- Y[-(1:2), ] - fitted1lev
  res1diff <- Y_diff_temp - fitted1diff
  res1 <- rbind(res1diff, res1lev)
  Z2 <- matrix(0, nrow = ncol(Z1), ncol = ncol(Z1))
  for (i in 1:N) {
    Z2 <- Z2 + t(as.matrix(S_i[[i]])) %*% (res1[, i] %*% t(res1[, i])) %*% S_i[[i]]
  }
  Z2try <- try(Z2 <- solve(Z2))
  if ("try-error" %in% class(Z2try)) {
    print("Second-step matrix is singular. Generalised inverse was used.")
    Z2 <- as.matrix(Z2)
    Z2 <- ginv(Z2)
  }
  rho2 <- solve(t(y_stacked_lag) %*% S %*% Z2 %*% t(as.matrix(S)) %*% y_stacked_lag) %*%
    (t(y_stacked_lag) %*% S %*% Z2 %*% t(as.matrix(S)) %*% y_stacked)
  if (model == "twosteps") {
    return(as.numeric(rho2))
  }

  # step 3
  fitted2lev <- matrix(NA, nrow = T - 2, ncol = N)
  fitted2diff <- matrix(NA, nrow = T - 2, ncol = N)
  for (t in 1:(T - 2)) {
    fitted2lev[t, ] <- Y[t + 1, ] * as.numeric(rho2)
    fitted2diff[t, ] <- Y_diff_temp[t, ] * as.numeric(rho2)
  }
  res2lev <- Y[-(1:2), ] - fitted2lev
  res2diff <- Y_diff_temp - fitted2diff
  res2 <- rbind(res2diff, res2lev)
  Z3 <- matrix(0, nrow = ncol(Z1), ncol = ncol(Z1))
  for (i in 1:N) {
    Z3 <- Z3 + t(as.matrix(S_i[[i]])) %*% (res2[, i] %*% t(res2[, i])) %*% S_i[[i]]
  }
  Z3try <- try(Z3 <- solve(Z3))
  if ("try-error" %in% class(Z3try)) {
    print("Third-step matrix is singular. Generalised inverse was used.")
    Z3 <- as.matrix(Z3)
    Z3 <- ginv(Z3)
  }
  rho3 <- solve(t(y_stacked_lag) %*% S %*% Z3 %*% t(as.matrix(S)) %*% y_stacked_lag) %*%
    (t(y_stacked_lag) %*% S %*% Z3 %*% t(as.matrix(S)) %*% y_stacked)
  if (model == "threesteps") {
    return(as.numeric(rho3))
  }
}
