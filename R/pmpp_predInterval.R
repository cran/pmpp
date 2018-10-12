#' Random-Window Block Bootstrap for prediction intervals for PMPP model
#'
#' @author Michal Oleszak
#'
#' @description Produces prediction intervals for Posterior Mean Panel Predictor
#'              model by means of resampling with replacement from model's residuals.
#'              Block Bootstrap method takes into account heteroskedasticity of the
#'              error terms, both across units and over time. Block window is chosen randomly.
#'
#' @references Oleszak, M. (2018). "Forecasting sales with micro-panels:
#'             Empirical Bayes approach. Evidence from consumer goods sector.",
#'             Erasmus University Thesis Repository
#'
#' @param model        PMPP model, as returned by \code{pmpp()}
#' @param boot_reps    integer; number of bootstrap replications
#' @param block_size   integer; width of the re-sampled block of residuals
#' @param confidence   numeric in (0,1); confidence level of the interval
#' @param iter         iterating constant, to be used in a loop when extraction
#'                     from call is needed
#' @param fframe       \code{data.frame} with the same columns as input data
#'                     to \code{model}, but with empty rows added to each
#'                     cross-sectional unit, as created by \code{create_fframe()}
#'
#' @return A \code{data.frame} with panel indices, lower and upper bounds and midpoint.
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom dplyr filter select bind_cols
#' @importFrom magrittr "%>%"
#' @export
#'
#' @examples
#' \dontrun{data(EmplUK, package = "plm")
#' EmplUK <- dplyr::filter(EmplUK, year %in% c(1978, 1979, 1980, 1981, 1982))
#' pmpp_model <- pmpp(dep_var = "emp", data = EmplUK)
#' my_fframe <- create_fframe(EmplUK, 1983:1985)
#' intervals <- pmpp_predinterval(pmpp_model, my_fframe, boot_reps = 10)
#' }
pmpp_predinterval <- function(model, fframe, boot_reps = 1000, block_size = NULL,
                              confidence = 0.95, iter = NULL) {
  start_global <- Sys.time()
  if (!is.null(iter)) {
    data_orig <- model$data[[iter]]
  } else {
    data_orig <- model$data
  }
  unit_id <- model$args$panel_ind[1]
  time_id <- model$args$panel_ind[2]
  old_time <- unique(data_orig[[time_id]])
  n.ahead <- sum(!(as.character(unique(fframe[[time_id]])) %in% old_time))
  N <- model$sample_size[1]
  T <- model$sample_size[2] - 1
  Y_all <- matrix(data_orig[[model$args$dep_var]], nrow = T + 1, ncol = N)
  if (length(model$common_coeff[-1]) > 0) {
    Z_all <- vector("list", length(model$common_coeff[-1]))
    for (z in seq(model$common_coeff[-1])) {
      Z_all[[z]] <- matrix(data_orig[[model$args$exp_var[z]]], nrow = T + 1, ncol = N)
    }
  }
  y_T <- Y_all[T + 1, ]

  # Specify default block size ------------------------------------------------
  if (is.null(block_size)) {
    block_size <- as.numeric(round(T ^ (1 / 3)))
  }

  # Extract data from model ---------------------------------------------------
  # (model's residuals u_hat and parameter estimates rho_hat, lambda_hat, alphas_hat)
  u_hat <- t(model$residuals)
  rho_hat <- model$common_coeff[1]
  lambda_hat <- model$intercept
  alphas_hat <- model$common_coeff[-1]

  # Initialise the bootstrap --------------------------------------------------
  bootstraps <- list()
  pb <- txtProgressBar(min = 0, max = boot_reps, style = 3, char = "~")
  for (R in 1:boot_reps) {

    # Create random blocks & stack them up ------------------------------------
    B <- matrix(nrow = N)
    while (ncol(B) - 1 < T) {
      j <- sample(1:(T - block_size + 1), 1, replace = TRUE)
      B_j <- u_hat[, j:(j + block_size - 1)]
      B <- cbind(B, B_j)
    }
    B <- B[, -1]
    if (ncol(B) > T) {
      B <- B[, 1:T]
    }

    # Generate bootstrap replicates y_boot ------------------------------------
    y_boot <- matrix(nrow = N, ncol = T)
    if (length(alphas_hat) > 0) {
      for (z in seq_along(Z_all)) {
        Z_all[[z]] <- alphas_hat[z] * Z_all[[z]]
      }
      Z_impact <- Reduce("+", Z_all)
    }
    for (t in 1:T) {
      if (length(alphas_hat) == 0) {
        y_boot[, t] <- lambda_hat + rho_hat * Y_all[t, ] + B[, t]
      } else {
        y_boot[, t] <- lambda_hat + rho_hat * Y_all[t, ] + Z_impact[t, ] + B[, t]
      }
    }

    # Re-fit the model to obtain bootstrap coefficients -----------------------
    # (lambda_boot, rho_boot, alpha_boot)
    y_boot <- pmpp_data(y_boot)
    if (!is.null(iter)) {
      indata <- model$data[[iter]]
    } else {
      indata <- model$data
    }
    new_index <- subset(indata[, 1:2], indata[2] != as.character(indata[1, 2]))
    newdata <- cbind(new_index, as.data.frame(y_boot)[3])
    newdata[, 1] <- as.character(newdata[, 1])
    newdata[, 3] <- as.numeric(newdata[, 3])
    if (length(alphas_hat) > 0) {
      exp_vars <- dplyr::filter(indata, get(time_id) %in% unique(newdata[[time_id]])) %>%
        dplyr::select(model$args$exp_var)
      newdata <- bind_cols(newdata, exp_vars)
    }
    bootmodel <- pmpp(
      dep_var = colnames(as.data.frame(y_boot))[3],
      panel_ind = colnames(new_index),
      exp_var = model$args$exp_var,
      csi_var = model$args$csi_var,
      data = newdata,
      post_mean_method = model$args$post_mean_method,
      common_par_method = model$args$common_par_method,
      optim_method = model$args$optim_method,
      dens_grid = model$args$dens_grid,
      gmm_model = model$args$gmm_model,
      gmm_inst = model$args$gmm_inst
    )
    rho_boot <- bootmodel$common_coeff[1]
    lambda_boot <- bootmodel$intercept
    alphas_boot <- bootmodel$common_coeff[-1]

    # Create random blocks for prediction & stack them up ---------------------
    B <- matrix(nrow = N)
    while (ncol(B) - 1 < n.ahead) {
      j <- sample(1:(T - block_size + 1), 1, replace = TRUE)
      B_j <- u_hat[, j:(j + block_size - 1)]
      B <- cbind(B, B_j)
    }
    B <- as.matrix(B[, -1])
    if (ncol(B) > n.ahead) {
      B <- B[, 1:n.ahead]
    }

    # Compute forecasts -------------------------------------------------------
    y_pred <- matrix(nrow = N, ncol = n.ahead)
    rho_vec <- c()
    for (i in 0:(n.ahead - 1)) {
      rho_vec[i + 1] <- rho_boot ^ i
    }
    # AR(1) model, no external variables
    if (length(alphas_boot) == 0) {
      for (h in 1:n.ahead) {
        y_pred[, h] <- lambda_boot * sum(rho_vec[1:h]) + rho_boot ^ h * y_T + B[, h]
      }
    } else {
      # ARX(1) model with predicion horizon = T+1
      if (n.ahead == 1) {
        Z_boot_impact <- vector("list", length(alphas_boot))
        names(Z_boot_impact) <- names(alphas_boot)
        data_last <- dplyr::filter(data_orig, get(time_id) == max(get(time_id)))
        for (a in names(alphas_boot)) {
          Z_boot_impact[[a]] <- alphas_boot[a] * data_last[[a]]
        }
        Z_boot_impact_all_exp_vars <- Reduce("+", Z_boot_impact)
        y_pred[, 1] <- lambda_boot + (rho_boot * y_T) + Z_boot_impact_all_exp_vars + B[, 1]
        # ARX(1) model with horizon > T+1
      } else {
        laststamp <- as.character(rev(unique(fframe[[time_id]]))[1])
        if (any(is.na(
          dplyr::filter(fframe, !get(time_id) %in% c(old_time, laststamp)) %>%
            select(model$args$exp_var)
        ))) {
          stop("Forecasting allowed only one step ahead further than the latest exp_vars. Adjust fframe.")
        }
        newdata1 <- dplyr::filter(fframe, get(time_id) != laststamp)
        new_time <- setdiff(unique(newdata1[[time_id]]), old_time)
        newdata2 <- dplyr::filter(newdata1, get(time_id) %in%
          c(new_time, max(data_orig[[time_id]])))
        T_and_later <- as.character(unique(newdata2[[time_id]]))
        Z_boot_impact <- vector("list", length(alphas_boot))
        names(Z_boot_impact) <- names(alphas_boot)
        Z_boot_impact <- lapply(Z_boot_impact, function(x) {
          x <- matrix(NA, ncol = n.ahead, nrow = length(lambda_boot))
        })
        for (a in names(alphas_boot)) {
          for (h in seq(n.ahead)) {
            Z_list <- vector("list", h)
            for (g in 1:h) {
              Z_list[[g]] <- dplyr::filter(newdata2, get(time_id) == T_and_later[g]) %>%
                dplyr::select(a)
            }
            for (z in seq(Z_list)) {
              Z_list[[z]] <- Z_list[[z]] * rev(rho_vec[1:h])[z]
            }
            sum_h <- Reduce("+", Z_list)
            Z_boot_impact[[a]][, h] <- as.matrix(alphas_boot[a] * sum_h)
          }
        }
        Z_boot_impact_all_exp_vars <- Reduce("+", Z_boot_impact)
        for (h in 1:n.ahead) {
          y_pred[, h] <- lambda_boot * sum(rho_vec[1:h]) + rho_boot ^ h * y_T +
            Z_boot_impact_all_exp_vars[, h] + B[, h]
        }
      }
    }

    # Save replication's outcome ----------------------------------------------
    bootstraps[[R]] <- y_pred
    setTxtProgressBar(pb, R)
  }
  close(pb)

  # Obtain quantiles ----------------------------------------------------------
  bootstraps <- lapply(bootstraps, pmpp_data)
  out <- matrix(nrow = N * n.ahead, ncol = 5)
  colnames(out) <- c(
    unit_id, time_id,
    paste0("lower_", 100 * confidence),
    paste0("upper_", 100 * confidence), "midpoint"
  )
  out <- as.data.frame(out)
  out[, 1:2] <- bootstraps[[1]][, 1:2]
  cutoff <- ceiling((boot_reps * (1 - confidence)) / 2)
  for (i in 1:(N * n.ahead)) {
    temp <- c()
    for (R in 1:boot_reps) {
      temp <- c(temp, bootstraps[[R]][i, 3])
    }
    temp <- c(sort(temp))
    temp <- temp[-c((length(temp):(length(temp) - cutoff + 1)), 1:cutoff)]
    out[[paste0("lower_", 100 * confidence)]][i] <- min(temp)
    out[[paste0("upper_", 100 * confidence)]][i] <- max(temp)
    out[["midpoint"]][i] <- mean(temp)
  }
  out[[unit_id]] <- factor(out[[unit_id]],
    levels = unique(model$data[[unit_id]])
  )
  levels(out[[time_id]]) <- setdiff(unique(fframe[[time_id]]), old_time)

  # Return the output ---------------------------------------------------------
  message(paste0(
    "Bootstrapping done! It took ",
    substring(capture.output(Sys.time() - start_global), 20), "."
  ))
  return(out)
}
