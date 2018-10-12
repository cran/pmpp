#' Posterior Mean Panel Predictor for dynamic panel modelling
#'
#' @author Michal Oleszak
#'
#' @description This function estimates parameters of the Posterior Mean Panel
#'              Predictor (PMPP) model based on an empirical-Bayes approach to
#'              obtain unit-specific fixed effects.
#'
#' @references Liu et al. (2016), "Forecasting with Dynamic Panel Data Models",
#'             PIER Working Paper No. 16022.,
#'             \url{https://papers.ssrn.com/sol3/Papers.cfm?abstract_id=2889000}
#' @references Oleszak, M. (2018). "Forecasting sales with micro-panels:
#'             Empirical Bayes approach. Evidence from consumer goods sector.",
#'             Erasmus University Thesis Repository
#'
#' @param dep_var character string indicating name of dependent variable
#' @param panel_ind vector of length 2 indicating names of variables indexing units and time periods respectively
#' @param exp_var vector of character strings indicating names of exogeneous explanatory variables
#' @param csi_var vector of character strings indicating names of cross-sectionally invariant explanatory variables; 
#' feature not supported yet
#' @param data data.frame or matrix with input data
#' @param post_mean_method method for estimating the heterogeneous intercept parameters, one of "gaussian", "kernel"
#' @param common_par_method method for estimating the common parameters, one of "QMLE", "GMM_ABond", "GMM_BBond", 
#' GMM_ABover", "GMM_SSYS"
#' @param optim_method which optimisation routine to use, one of "gradient", "quadratic", "annealing"
#' @param dens_grid size of the grid over which data is interpolated for kernel density estimation; 
#' larger value may yield higher accuracy, but increases computation time
#' @param gmm_model number of steps for computing optimal GMM matrix, one of "onestep", "twosteps", 
#' "threesteps"; "threesteps" can be used for "GMM_SSYS" only
#' @param gmm_inst number of lagged values of the dependent variable to be used as GMM instruments 
#' in Arellano-Bond/Blundell-Bond setting
#' @param pure_data if TRUE, removes indexing/subsetting from model's call on data, facilitating use in a loop
#'
#' @details The PMPP model is a two-step procedure. First, the homogeneous parameters are
#'          estimated using one of the QMLE or GMM-based methods:
#'          \itemize{
#'          \item Arellano-Bond estimator (Difference GMM),
#'          \item Arellano-Bover estimator (Level GMM),
#'          \item Blundell-Bond estimator (System GMM),
#'          \item Sub-optimal System GMM estimator,
#'          \item Quasi-Maximum Likelihood estimator.
#'          }
#'          Parameter \code{common_par_method} can be used to select the method for common parameters estimation.
#'          All the above methods only provide estimates of the homogeneous parameters, i.e. the
#'          ones measuring impact of lagged response and external variables. The intercept is removed in the
#'          estimation process.
#'          In the second step of the PMPP modelling, the individual-specific intercept is
#'          calculated based on the formula for posterior mean (Tweedie's Formula). It involves
#'          approximating certain density function, which can be done in two ways:
#'          \itemize{
#'          \item Parametrically, assuming Gaussian distribution,
#'          \item Using a 2D kernel density estimator.
#'          }
#'          Parameter \code{post_mean_method} can be used to select the method used for intercept estimation.
#'          For technical details on the methods, see the references.
#'
#' @return An object of class \code{pmpp}; a list with parameter estimates, fitted values,
#'         residuals, in-sample error measures and information on the data and function call.
#'
#' @importFrom pracma mrdivide mldivide zeros ones
#' @importFrom plm plm pgmm
#' @importFrom minqa bobyqa
#' @importFrom data.table setattr
#' @importFrom utils capture.output tail
#' @importFrom stats optim as.formula var
#' @export
#'
#' @examples
#' data(EmplUK, package = "plm")
#' EmplUK <- dplyr::filter(EmplUK, year %in% c(1978, 1979, 1980, 1981, 1982))
#' pmpp_model <- pmpp(dep_var = "emp", data = EmplUK)
pmpp <- function(dep_var,
                 data,
                 panel_ind = colnames(data[, 1:2]),
                 exp_var = NULL,
                 csi_var = NULL,
                 post_mean_method = "gaussian",
                 common_par_method = "QMLE",
                 optim_method = "quadratic",
                 dens_grid = 2 ^ 10,
                 gmm_model = "twosteps",
                 gmm_inst = 99,
                 pure_data = FALSE) {
  start_global <- Sys.time()

  # Validate input ---------------------------------------------------------------
  # dep_var
  if (missing(dep_var) | class(dep_var) != "character") {
    stop("Dependent variable's name must be specified as a character string.")
  }
  # panel_ind
  if (class(panel_ind) != "character" | length(panel_ind) != 2) {
    stop("panel_ind must be specified as vector of length 2 with the first and
        second element being names of variables indexing units and time
        periods respectively.")
  }
  # exp_var
  if (class(exp_var) != "character" & !is.null(exp_var)) {
    stop("Explanatory variables' names must be specified as a vector of character
       strings.")
  }
  # csi_var
  if (class(csi_var) != "character" & !is.null(csi_var)) {
    stop("Cross-sectionally invariant variables' names must be specified as
        a vector of character strings.")
  }
  if (!is.null(csi_var)) {
    stop("This feature is not supported yet. Use default value: csi_var = NULL.")
  }
  # data
  if (missing(data) | !class(data) %in% c("data.frame", "matrix")) {
    "data must be of class 'data.frame' or 'matrix'"
  }
  # post_mean_method
  if (!post_mean_method %in% c("gaussian", "kernel")) {
    stop("post_mean_method must be one of: 'gaussian' or 'kernel'")
  }
  # common_par_method
  if (!common_par_method %in% c(
    "QMLE", "GMM_ABond", "GMM_BBond", "GMM_ABover",
    "GMM_SSYS"
  )) {
    stop("common_par_method must be one of: 'QMLE', 'GMM_ABond', 'GMM_BBond',
       'GMM_ABover', 'GMM_SSYS'")
  }
  # optim_method
  if (!optim_method %in% c("gradient", "quadratic", "annealing")) {
    stop("optim_method must be one of: 'gradient', 'quadratic', 'annealing'")
  }
  # gmm_model
  if (!gmm_model %in% c("onestep", "twosteps", "threesteps") |
    (gmm_model == "threesteps" & common_par_method != "GMM_SSYS")) {
    stop("gmm_model must be one of: 'onestep', 'twosteps', 'threesteps'.
        'threesteps' can be used for GMM_SSYS only.")
  }
  # gmm_inst
  if (!class(gmm_inst) %in% c("numeric", "integer") | gmm_inst < 1 |
    gmm_inst != round(gmm_inst)) {
    stop("gmm_inst must be a positive integer.")
  }
  # dens_grid
  if (!class(dens_grid) %in% c("numeric", "integer") | dens_grid < 1 |
    dens_grid != round(dens_grid)) {
    stop("dens_grid must be a positive integer.")
  }
  # pure_data
  if (class(pure_data) != "logical") {
    stop("pure_data must be logical")
  }

  # Data prep --------------------------------------------------------------------
  unit <- as.matrix(data[which(colnames(data) == panel_ind[1])])
  time <- as.matrix(data[which(colnames(data) == panel_ind[2])])
  dep_var <- as.matrix(data[which(colnames(data) == dep_var)])
  N <- as.numeric(nrow(unique(unit)))
  T_orig <- nrow(unique(time))
  if (is.null(exp_var)) {
    n_alpha <- 0
    indata <- data.frame(unit, time, dep_var)
  } else {
    exp_var <- as.matrix(data[which(colnames(data) %in% exp_var)])
    n_alpha <- ncol(exp_var)
    zname <- colnames(exp_var)
    indata <- data.frame(unit, time, dep_var, exp_var)
  }
  yname <- names(as.data.frame(dep_var))
  Y0 <- matrix(dep_var, nrow = T_orig, ncol = N)[1, ]
  Y_all <- matrix(dep_var, nrow = T_orig, ncol = N)
  Y_mat <- matrix(dep_var, nrow = T_orig, ncol = N)[2:T_orig, ]
  X_mat <- matrix(dep_var, nrow = T_orig, ncol = N)[1:T_orig - 1, ]
  y_T <- c(tail(Y_mat, 1))
  if (n_alpha > 0) {
    Z_mat <- vector("list", n_alpha)
    for (a in 1:n_alpha) {
      Z_mat[[a]] <- matrix(exp_var[, a], nrow = T_orig, ncol = N)[2:T_orig, ]
    }
  } else {
    Z_mat <- NULL
  }
  T <- nrow(Y_mat)
  if (is.null(csi_var)) {
    W <- cbind(rep(1, T))
  } else {
    W_mat <- matrix(data[[csi_var]], nrow = T_orig, ncol = N)
    W_check <- apply(W_mat, 1, function(x) length(unique(x)))
    if (any(W_check != 1)) {
      stop(paste0(csi_var, " is not cross-sectionally invariant!"))
    }
    W <- cbind(rep(1, T), W_mat[-1, 1])
  }
  n_lambda <- ncol(W)
  n_param <- 2 + (n_lambda * 3) + n_alpha
  # Orthogonal Forward Deviation transformation for Arellano-Bover estimator
  if (common_par_method == "GMM_ABover") {
    Y_star <- matrix(NA, nrow = T - n_lambda, ncol = N)
    X_star <- matrix(NA, nrow = T - n_lambda, ncol = N)
    if (!is.null(Z_mat)) {
      Z_star <- vector("list", n_alpha)
      for (a in 1:n_alpha) {
        Z_star[[a]] <- matrix(NA, nrow = T - n_lambda, ncol = N)
      }
    } else {
      Z_star <- NULL
    }
    for (t in 1:(T - n_lambda)) {
      c <- sqrt(1 + mrdivide(W[t, ], matrix(W[(t + 1):nrow(W), ], ncol = ncol(W))) %*%
        W[(t + 1):nrow(W), ] %*% W[t, ])
      Y_star[t, ] <- (Y_mat[t, ] - (mrdivide(W[t, ], matrix(W[(t + 1):nrow(W), ],
        ncol = ncol(W)
      )) %*% W[(t + 1):nrow(W), ]) %*% t(W[(t + 1):nrow(W), ]) %*%
        Y_mat[(t + 1):nrow(Y_mat), ]) / as.numeric(c)
      X_star[t, ] <- (X_mat[t, ] - (mrdivide(W[t, ], matrix(W[(t + 1):nrow(W), ],
        ncol = ncol(W)
      )) %*% W[(t + 1):nrow(W), ]) %*% t(W[(t + 1):nrow(W), ]) %*%
        X_mat[(t + 1):nrow(X_mat), ]) / as.numeric(c)
      if (!is.null(Z_star)) {
        for (a in 1:n_alpha) {
          Z_star[[a]][t, ] <- (Z_mat[[a]][t, ] - (mrdivide(
            W[t, ],
            matrix(W[(t + 1):nrow(W), ], ncol = ncol(W))
          ) %*%
            W[(t + 1):nrow(W), ]) %*% t(W[(t + 1):nrow(W), ]) %*%
            Z_mat[[a]][(t + 1):nrow(Z_mat[[a]]), ]) / as.numeric(c)
        }
      }
    }
  } else {
    X_star <- Y_star <- Z_star <- NULL
  }

  # Compute initial estimator ----------------------------------------------------
  if (n_alpha == 0) {
    dX <- (diag(T) - mrdivide(W, t(W) %*% W) %*% t(W)) %*% X_mat
    rho_0 <- mldivide(t(c(dX)) %*% c(dX), t(c(dX))) %*% c(Y_mat)
    lambda_0 <- get_lambda0(rho_0,
      N = N, T = T, n_alpha = n_alpha,
      Y_mat = Y_mat, X_mat = X_mat, W = W, Z_mat = Z_mat
    )
    sigma2_0 <- get_sigma2(rho_0,
      common_par_method = common_par_method, X_star = X_star,
      Y_star = Y_star, Z_star = Z_star, X_mat = X_mat,
      Y_mat = Y_mat, Z_mat = Z_mat, n_alpha = n_alpha
    )
  } else {
    fmla <- as.formula(paste(
      paste(yname, " ~ ", collapse = " "),
      paste0("lag(", yname, ",1) + "), paste(zname, collapse = " + ")
    ))
    model_init <- plm(fmla, data = indata, model = "within")
    rho_0 <- model_init$coefficients[1]
    alpha_0 <- c(model_init$coefficients[2:(n_alpha + 1)])
    lambda_0 <- get_lambda0(rho_0, alpha_0,
      N = N, T = T, n_alpha = n_alpha,
      Y_mat = Y_mat, X_mat = X_mat, W = W, Z_mat = Z_mat
    )
    sigma2_0 <- get_sigma2(rho_0, alpha_0,
      common_par_method = common_par_method,
      X_star = X_star, Y_star = Y_star, Z_star = Z_star,
      X_mat = X_mat, Y_mat = Y_mat, Z_mat = Z_mat,
      n_alpha = n_alpha
    )
  }

  # Initialise optimisation routine ----------------------------------------------
  init <- zeros(n_param, 1)
  if (n_alpha == 0) {
    init[1:2] <- c(rho_0, sigma2_0)
  } else {
    init[1:(n_alpha + 2)] <- c(rho_0, alpha_0, sigma2_0)
  }
  aux_Y0 <- rbind(ones(1, N), Y0)
  aux_coef <- mldivide(aux_Y0 %*% t(aux_Y0), aux_Y0) %*% t(lambda_0)
  init[(n_alpha + 3):(length(init) - n_lambda)] <- t(c(aux_coef))
  aux_res <- t(lambda_0) - t(aux_Y0) %*% aux_coef
  init[(length(init) - n_lambda + 1):length(init)] <- diag(var(aux_res))

  # QMLE -------------------------------------------------------------------------
  if (common_par_method == "QMLE") {
    if (optim_method == "gradient") {
      minimise_QMLE <- try(
        minimise_QMLE <- optim(init, loglikelihood_QMLE,
          n_alpha = n_alpha,
          X_mat = X_mat, Y_mat = Y_mat,
          Z_mat = Z_mat, W = W, N = N,
          T = T, aux_Y0 = aux_Y0
        )
      )
      if ("try-error" %in% class(minimise_QMLE)) {
        print("Gradient optimisation failed. Quadratic approximation was used.")
        minimise_QMLE <- bobyqa(init, loglikelihood_QMLE,
          n_alpha = n_alpha,
          X_mat = X_mat, Y_mat = Y_mat, Z_mat = Z_mat,
          W = W, N = N, T = T, aux_Y0 = aux_Y0
        )
      }
    } else if (optim_method == "quadratic") {
      minimise_QMLE <- bobyqa(init, loglikelihood_QMLE,
        n_alpha = n_alpha,
        X_mat = X_mat, Y_mat = Y_mat, Z_mat = Z_mat,
        W = W, N = N, T = T, aux_Y0 = aux_Y0
      )
    } else if (optim_method == "annealing") {
      minimise_QMLE <- optim(init, loglikelihood_QMLE,
        n_alpha = n_alpha,
        X_mat = X_mat, Y_mat = Y_mat, Z_mat = Z_mat, W = W,
        N = N, T = T, aux_Y0 = aux_Y0, method = "SANN"
      )
    }
    param_QMLE <- minimise_QMLE$par
    rho_QMLE_par <- param_QMLE[1]
    if (n_alpha > 0) {
      alpha_QMLE_par <- param_QMLE[2:(n_alpha + 1)]
    }
    sigma2_QMLEpar <- param_QMLE[n_alpha + 2]
    mmu_QMLEpar <- matrix(param_QMLE[(n_alpha + 3):
    (length(param_QMLE) - n_lambda)], nrow = 2)
    ww2_lambda_QMLEpar <- param_QMLE[(length(param_QMLE) - n_lambda + 1):
    length(param_QMLE)]
    lambda_QMLEpar <- get_lambda0(rho_QMLE_par,
      N = N, T = T,
      n_alpha = n_alpha, Y_mat = Y_mat,
      X_mat = X_mat, W = W, Z_mat = Z_mat
    )
    if (post_mean_method == "gaussian") {
      post_mean_method_lambda_QMLE_par <- post_mean_lambda_par(lambda_QMLEpar,
        sigma2_QMLEpar, mmu_QMLEpar, ww2_lambda_QMLEpar,
        W = W, aux_Y0 = aux_Y0
      )
    } else if (post_mean_method == "kernel") {
      post_mean_method_lambda_QMLE_kernel <- suppressWarnings(
        get_kernel(lambda_QMLEpar, sigma2_QMLEpar,
          dens_grid = dens_grid,
          N = N, T = T, Y0 = Y0
        )
      )
    }
  }

  # GMM Arellano&Bond ------------------------------------------------------------
  if (common_par_method == "GMM_ABond") {
    if (n_alpha == 0) {
      fmla <- as.formula(paste(
        paste(yname, " ~ ", collapse = " "),
        paste0("lag(", yname, ", 1)"),
        paste0(" | lag(", paste0(yname, ", 2:", gmm_inst, ")"))
      ))
      gmmab <- pgmm(fmla,
        data = indata,
        effect = "individual", model = gmm_model, transformation = "d"
      )
      if (gmm_model == "onestep") {
        rho_ABond_par <- matrix(gmmab$coefficients)
      } else {
        rho_ABond_par <- matrix(gmmab$coefficients[[2]])
      }
      rho_GMMpar <- rho_ABond_par
      sigma2_GMMpar <- c(get_sigma2(rho_ABond_par,
        common_par_method = common_par_method,
        X_star = X_star, Y_star = Y_star,
        Z_star = Z_star, X_mat = X_mat, Y_mat = Y_mat,
        Z_mat = Z_mat, n_alpha = n_alpha
      ))
      if (post_mean_method == "gaussian") {
        post_mean_method_lambda_ABond_par <- GMM_parametric(rho_ABond_par,
          optim_method = optim_method,
          init = init,
          n_lambda = n_lambda,
          n_alpha = n_alpha,
          X_mat = X_mat, Y_mat = Y_mat,
          Z_mat = Z_mat, W = W, T = T,
          N = N, aux_Y0 = aux_Y0,
          common_par_method = common_par_method,
          X_star = X_star,
          Y_star = Y_star,
          Z_star = Z_star
        )
      } else if (post_mean_method == "kernel") {
        lambda_ABondpar <- get_lambda0(rho_ABond_par,
          N = N, T = T,
          n_alpha = n_alpha, Y_mat = Y_mat,
          X_mat = X_mat, W = W, Z_mat = Z_mat
        )
        sigma2_ABondpar <- c(get_sigma2(rho_ABond_par,
          common_par_method = common_par_method,
          X_star = X_star, Y_star = Y_star,
          Z_star = Z_star, X_mat = X_mat,
          Y_mat = Y_mat, Z_mat = Z_mat,
          n_alpha = n_alpha
        ))
        post_mean_method_lambda_ABond_kernel <- suppressWarnings(
          get_kernel(lambda_ABondpar,
            sigma2_ABondpar,
            dens_grid = dens_grid,
            N = N, T = T, Y0 = Y0
          )
        )
      }
    } else {
      fmla <- as.formula(paste(
        paste(yname, " ~ ", collapse = " "),
        paste0("lag(", yname, ", 1)", " + "),
        paste(zname, collapse = " + "),
        paste0(" | lag(", paste0(yname, ", 2:", gmm_inst, ")"))
      ))
      gmmab <- pgmm(fmla,
        data = indata,
        effect = "individual", model = gmm_model, transformation = "d"
      )
      if (gmm_model == "onestep") {
        rho_ABond_par <- matrix(gmmab$coefficients)[1]
        alpha_ABond_par <- matrix(gmmab$coefficients)[2:(n_alpha + 1)]
      } else {
        rho_ABond_par <- matrix(gmmab$coefficients[[2]])[1]
        alpha_ABond_par <- matrix(gmmab$coefficients[[2]])[2:(n_alpha + 1)]
      }
      rho_GMMpar <- rho_ABond_par
      alpha_GMMpar <- alpha_ABond_par
      sigma2_GMMpar <- c(get_sigma2(rho_ABond_par,
        alpha = alpha_ABond_par,
        common_par_method = common_par_method, X_star = X_star,
        Y_star = Y_star, Z_star = Z_star,
        X_mat = X_mat, Y_mat = Y_mat, Z_mat = Z_mat,
        n_alpha = n_alpha
      ))
      if (post_mean_method == "gaussian") {
        post_mean_method_lambda_ABond_par <- GMM_parametric(rho_ABond_par,
          alpha_ABond_par,
          optim_method = optim_method,
          init = init,
          n_lambda = n_lambda,
          n_alpha = n_alpha,
          X_mat = X_mat, Y_mat = Y_mat,
          Z_mat = Z_mat, W = W, T = T,
          N = N, aux_Y0 = aux_Y0,
          common_par_method = common_par_method,
          X_star = X_star,
          Y_star = Y_star,
          Z_star = Z_star
        )
      } else if (post_mean_method == "kernel") {
        lambda_ABondpar <- get_lambda0(rho_ABond_par,
          N = N, T = T,
          n_alpha = n_alpha, Y_mat = Y_mat,
          X_mat = X_mat, W = W, Z_mat = Z_mat
        )
        sigma2_ABondpar <- c(get_sigma2(rho_ABond_par,
          alpha = alpha_ABond_par,
          common_par_method = common_par_method, X_star = X_star,
          Y_star = Y_star, Z_star = Z_star,
          X_mat = X_mat, Y_mat = Y_mat,
          Z_mat = Z_mat, n_alpha = n_alpha
        ))
        post_mean_method_lambda_ABond_kernel <- suppressWarnings(
          get_kernel(lambda_ABondpar,
            sigma2_ABondpar,
            dens_grid = dens_grid,
            N = N, T = T, Y0 = Y0
          )
        )
      }
    }
  }

  # GMM Blundell&Bond ------------------------------------------------------------
  if (common_par_method == "GMM_BBond") {
    if (n_alpha == 0) {
      fmla <- as.formula(paste(
        paste(yname, " ~ ", collapse = " "),
        paste0("lag(", yname, ", 1)"),
        paste0(" | lag(", paste0(yname, ", 2:", gmm_inst, ")"))
      ))
      gmmbb <- pgmm(fmla,
        data = indata,
        effect = "individual", model = gmm_model, transformation = "ld"
      )
      if (gmm_model == "onestep") {
        rho_BBond_par <- matrix(gmmbb$coefficients)
      } else {
        rho_BBond_par <- matrix(gmmbb$coefficients[[2]])
      }
      rho_GMMpar <- rho_BBond_par
      sigma2_GMMpar <- c(get_sigma2(rho_BBond_par,
        common_par_method = common_par_method,
        X_star = X_star, Y_star = Y_star,
        Z_star = Z_star, X_mat = X_mat, Y_mat = Y_mat,
        Z_mat = Z_mat, n_alpha = n_alpha
      ))
      if (post_mean_method == "gaussian") {
        post_mean_method_lambda_BBond_par <- GMM_parametric(rho_BBond_par,
          optim_method = optim_method,
          init = init,
          n_lambda = n_lambda,
          n_alpha = n_alpha,
          X_mat = X_mat, Y_mat = Y_mat,
          Z_mat = Z_mat, W = W, T = T,
          N = N, aux_Y0 = aux_Y0,
          common_par_method = common_par_method,
          X_star = X_star,
          Y_star = Y_star,
          Z_star = Z_star
        )
      } else if (post_mean_method == "kernel") {
        lambda_BBondpar <- get_lambda0(rho_BBond_par,
          N = N, T = T,
          n_alpha = n_alpha, Y_mat = Y_mat,
          X_mat = X_mat, W = W, Z_mat = Z_mat
        )
        sigma2_BBondpar <- c(get_sigma2(rho_BBond_par,
          common_par_method = common_par_method,
          X_star = X_star, Y_star = Y_star,
          Z_star = Z_star, X_mat = X_mat,
          Y_mat = Y_mat, Z_mat = Z_mat,
          n_alpha = n_alpha
        ))
        post_mean_method_lambda_BBond_kernel <- suppressWarnings(
          get_kernel(lambda_BBondpar,
            sigma2_BBondpar,
            dens_grid = dens_grid,
            N = N, T = T, Y0 = Y0
          )
        )
      }
    } else {
      fmla <- as.formula(paste(
        paste(yname, " ~ ", collapse = " "),
        paste0("lag(", yname, ", 1)", " + "),
        paste(zname, collapse = " + "),
        paste0(" | lag(", paste0(yname, ", 2:", gmm_inst, ")"))
      ))
      gmmbb <- pgmm(fmla,
        data = indata,
        effect = "individual", model = gmm_model, transformation = "ld"
      )
      if (gmm_model == "onestep") {
        rho_BBond_par <- matrix(gmmbb$coefficients)[1]
        alpha_BBond_par <- matrix(gmmbb$coefficients)[2:(n_alpha + 1)]
      } else {
        rho_BBond_par <- matrix(gmmbb$coefficients[[2]])[1]
        alpha_BBond_par <- matrix(gmmbb$coefficients[[2]])[2:(n_alpha + 1)]
      }
      rho_GMMpar <- rho_BBond_par
      alpha_GMMpar <- alpha_BBond_par
      sigma2_GMMpar <- c(get_sigma2(rho_BBond_par,
        alpha = alpha_BBond_par,
        common_par_method = common_par_method, X_star = X_star,
        Y_star = Y_star, Z_star = Z_star,
        X_mat = X_mat, Y_mat = Y_mat, Z_mat = Z_mat,
        n_alpha = n_alpha
      ))
      if (post_mean_method == "gaussian") {
        post_mean_method_lambda_BBond_par <- GMM_parametric(rho_BBond_par,
          alpha = alpha_BBond_par,
          optim_method = optim_method,
          init = init,
          n_lambda = n_lambda,
          n_alpha = n_alpha,
          X_mat = X_mat, Y_mat = Y_mat,
          Z_mat = Z_mat, W = W, T = T,
          N = N, aux_Y0 = aux_Y0,
          common_par_method = common_par_method,
          X_star = X_star,
          Y_star = Y_star,
          Z_star = Z_star
        )
      } else if (post_mean_method == "kernel") {
        lambda_BBondpar <- get_lambda0(rho_BBond_par,
          N = N, T = T,
          n_alpha = n_alpha, Y_mat = Y_mat,
          X_mat = X_mat, W = W, Z_mat = Z_mat
        )
        sigma2_BBondpar <- c(get_sigma2(rho_BBond_par,
          alpha = alpha_BBond_par,
          common_par_method = common_par_method, X_star = X_star,
          Y_star = Y_star, Z_star = Z_star,
          X_mat = X_mat, Y_mat = Y_mat,
          Z_mat = Z_mat, n_alpha = n_alpha
        ))
        post_mean_method_lambda_BBond_kernel <- suppressWarnings(
          get_kernel(lambda_BBondpar,
            sigma2_BBondpar,
            dens_grid = dens_grid,
            N = N, T = T, Y0 = Y0
          )
        )
      }
    }
  }

  # GMM Arellano&Bover -----------------------------------------------------------
  if (common_par_method == "GMM_ABover") {
    aux_x <- 0
    aux_y <- 0
    if (!is.null(Z_star)) {
      aux_z <- list(0, length(Z_star))
    }
    for (t in 1:(T - n_lambda)) {
      L <- matrix(X_mat[1:t, ], ncol = ncol(X_mat))
      aux_x <- aux_x + mrdivide(X_star[t, ] %*% t(L), L %*% t(L)) %*% L %*% X_star[t, ]
      aux_y <- aux_y + mrdivide(X_star[t, ] %*% t(L), L %*% t(L)) %*% L %*% Y_star[t, ]
      if (!is.null(Z_star)) {
        for (a in 1:length(Z_star)) {
          aux_z[[a]] <- aux_z[[a]] + mrdivide(X_star[t, ] %*% t(L), L %*% t(L)) %*%
            L %*% Z_star[[a]][t, ]
        }
      }
    }
    if (n_alpha == 0) {
      rho_ABover_par <- mldivide(aux_x, aux_y)
      rho_GMMpar <- rho_ABover_par
      sigma2_GMMpar <- c(get_sigma2(rho_ABover_par,
        common_par_method = common_par_method,
        X_star = X_star, Y_star = Y_star,
        Z_star = Z_star, X_mat = X_mat, Y_mat = Y_mat,
        Z_mat = Z_mat, n_alpha = n_alpha
      ))
      if (post_mean_method == "gaussian") {
        post_mean_method_lambda_ABover_par <- GMM_parametric(rho_ABover_par,
          optim_method = optim_method,
          init = init,
          n_lambda = n_lambda,
          n_alpha = n_alpha,
          X_mat = X_mat, Y_mat = Y_mat,
          Z_mat = Z_mat, W = W, T = T,
          N = N, aux_Y0 = aux_Y0,
          common_par_method = common_par_method,
          X_star = X_star,
          Y_star = Y_star,
          Z_star = Z_star
        )
      } else if (post_mean_method == "kernel") {
        lambda_ABoverpar <- get_lambda0(rho_ABover_par,
          N = N, T = T,
          n_alpha = n_alpha, Y_mat = Y_mat,
          X_mat = X_mat, W = W, Z_mat = Z_mat
        )
        sigma2_ABoverpar <- c(get_sigma2(rho_ABover_par,
          common_par_method = common_par_method,
          X_star = X_star, Y_star = Y_star,
          Z_star = Z_star, X_mat = X_mat,
          Y_mat = Y_mat, Z_mat = Z_mat,
          n_alpha = n_alpha
        ))
        post_mean_method_lambda_ABover_kernel <- suppressWarnings(
          get_kernel(lambda_ABoverpar,
            sigma2_ABoverpar,
            dens_grid = dens_grid,
            N = N, T = T, Y0 = Y0
          )
        )
      }
    } else {
      stop("Arellano&Bover estimation is not supported in the presence of explanatory variables yet.")
    }
  }

  # SSYS.GMM --------------------------------------------------------------------
  if (common_par_method == "GMM_SSYS") {
    if (n_alpha == 0) {
      rho_SSYS_par <- ssys_gmm(Y_all, model = gmm_model)
      rho_GMMpar <- rho_SSYS_par
      sigma2_GMMpar <- c(get_sigma2(rho_SSYS_par,
        common_par_method = common_par_method,
        X_star = X_star, Y_star = Y_star,
        Z_star = Z_star, X_mat = X_mat,
        Y_mat = Y_mat, Z_mat = Z_mat,
        n_alpha = n_alpha
      ))
      if (post_mean_method == "gaussian") {
        post_mean_method_lambda_SSYS_par <- GMM_parametric(rho_SSYS_par,
          optim_method = optim_method,
          init = init,
          n_lambda = n_lambda,
          n_alpha = n_alpha,
          X_mat = X_mat, Y_mat = Y_mat,
          Z_mat = Z_mat, W = W, T = T,
          N = N, aux_Y0 = aux_Y0,
          common_par_method = common_par_method,
          X_star = X_star,
          Y_star = Y_star,
          Z_star = Z_star
        )
      } else if (post_mean_method == "kernel") {
        lambda_SSYSpar <- get_lambda0(rho_SSYS_par,
          N = N, T = T,
          n_alpha = n_alpha, Y_mat = Y_mat,
          X_mat = X_mat, W = W, Z_mat = Z_mat
        )
        sigma2_SSYSpar <- c(get_sigma2(rho_SSYS_par,
          common_par_method = common_par_method,
          X_star = X_star, Y_star = Y_star,
          Z_star = Z_star, X_mat = X_mat,
          Y_mat = Y_mat, Z_mat = Z_mat,
          n_alpha = n_alpha
        ))
        post_mean_method_lambda_SSYS_kernel <- suppressWarnings(
          get_kernel(lambda_SSYSpar,
            sigma2_SSYSpar,
            dens_grid = dens_grid,
            N = N, T = T, Y0 = Y0
          )
        )
      }
    } else {
      stop("SSYS GMM estimation is not supported in the presence of explanatory variables yet.")
    }
  }

  # Reutrn the output ------------------------------------------------------------
  # Return parameters
  intercept_options <- c(
    "post_mean_method_lambda_QMLE_par",
    "post_mean_method_lambda_QMLE_kernel",
    "post_mean_method_lambda_ABond_par",
    "post_mean_method_lambda_ABond_kernel",
    "post_mean_method_lambda_BBond_par",
    "post_mean_method_lambda_BBond_kernel",
    "post_mean_method_lambda_ABover_par",
    "post_mean_method_lambda_ABover_kernel",
    "post_mean_method_lambda_SSYS_par",
    "post_mean_method_lambda_SSYS_kernel"
  )
  rho_options <- c(
    "rho_QMLE_par", "rho_ABond_par", "rho_BBond_par",
    "rho_ABover_par", "rho_SSYS_par"
  )
  alpha_options <- c(
    "alpha_QMLE_par", "alpha_ABond_par", "alpha_BBond_par",
    "alpha_ABover_par"
  )
  intercept <- get(intercept_options[sapply(intercept_options, exists,
    envir = environment()
  )])
  rho <- get(rho_options[sapply(rho_options, exists,
    envir = environment()
  )])
  if (n_alpha > 0) {
    alphas <- get(alpha_options[sapply(alpha_options, exists,
      envir = environment()
    )])
    common_coeff <- c(rho, alphas)
    names(common_coeff) <- c(paste0("lag_", yname), zname)
  } else {
    common_coeff <- c(rho)
    names(common_coeff) <- c(paste0("lag_", yname))
  }
  # Return fitted values & residuals
  fitted_y <- matrix(ncol = N, nrow = T)
  Z_impact <- matrix(0, ncol = N, nrow = T)
  if (!is.null(Z_mat)) {
    for (a in 1:n_alpha) {
      Z_impact <- Z_impact + common_coeff[a + 1] * Z_mat[[a]]
    }
  }
  for (t in 1:T) {
    fitted_y[t, ] <- intercept + common_coeff[1] * Y_all[t, ] + Z_impact[t, ]
  }
  residuals <- Y_all[-1, ] - fitted_y
  # Return RMSE
  RMSE <- sqrt(sum(residuals ^ 2) / (N * T))
  # Return MAPE
  MAPE <- mean(abs(residuals / Y_all[-1, ]))
  # Return sample size
  sample_size <- c("N" = N, "T" = T_orig)
  # Return arguments
  args <- mget(names(formals()), sys.frame(sys.nframe()))
  args$data <- deparse(substitute(data))
  if (pure_data == TRUE & grepl("\\[\\[.\\]\\]", args$data)) {
    args$data <- substr(args$data, 1, regexpr("\\[\\[.\\]\\]", args$data)[[1]] - 1)
  }
  args$dep_var <- yname
  if (n_alpha > 0) {
    args$exp_var <- zname
  } else {
    args$exp_var <- NULL
  }
  # Yield output list
  output <- list(
    "intercept" = intercept,
    "common_coeff" = common_coeff,
    "fitted_values" = fitted_y,
    "residuals" = residuals,
    "sample_size" = sample_size,
    "RMSE" = RMSE,
    "MAPE" = MAPE,
    "data" = data,
    "args" = args,
    "call" = sys.calls()[[1]],
    "elapsed_time" = substring(capture.output(Sys.time() - start_global), 20)
  )
  data.table::setattr(output, "class", c("pmpp"))
  return(output)
}
