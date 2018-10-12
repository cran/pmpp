#' Compute forecasts with a PMPP model
#'
#' @author Michal Oleszak
#'
#' @param object   an object of class \code{pmpp()}
#' @param fframe   \code{data.frame} with the same columns as input data
#'                 to \code{model}, but with empty rows added to each
#'                 cross-sectional unit, as created by \code{create_fframe()}
#' @param iter     iterating constant, to be used in a loop when extraction
#'                 from call is needed
#' @param ...      other arguments passed to the method
#'
#' @return A \code{data.frame} with predicted and true values.
#'
#' @importFrom dplyr mutate select left_join filter setdiff
#' @importFrom magrittr "%>%"
#' @importFrom utils head
#' @export
#'
#' @examples
#' data(EmplUK, package = "plm")
#' EmplUK <- dplyr::filter(EmplUK, year %in% c(1978, 1979, 1980, 1981, 1982))
#' pmpp_model <- pmpp(dep_var = "emp", data = EmplUK)
#' my_fframe <- create_fframe(EmplUK, 1983:1985)
#' prediction <- predict(pmpp_model, my_fframe)
predict.pmpp <- function(object, fframe = NULL, iter = NULL, ...) {
  A <- B <- NULL
  model <- object
  if (!inherits(model, "pmpp")) {
    stop("Non-convenient object, model should be of class 'pmpp'.")
  }

  # If no fframe provided, return fitted values -------------------------------
  if (is.null(fframe)) {
    return(model$fitted_values)
  }

  # Get model output & prepare forecast frame ---------------------------------
  if (!is.null(iter)) {
    y_T <- matrix(model$data[[iter]][[model$args$dep_var]],
      nrow = model$sample_size[2],
      ncol = model$sample_size[1]
    )[model$sample_size[2], ]
    data_orig <- model$data[[iter]]
  } else {
    y_T <- matrix(model$data[[model$args$dep_var]],
      nrow = model$sample_size[2],
      ncol = model$sample_size[1]
    )[model$sample_size[2], ]
    data_orig <- model$data
  }
  time_id <- model$args$panel_ind[2]
  unit_id <- model$args$panel_ind[1]
  laststamp <- as.character(rev(unique(fframe[[time_id]]))[1])
  old_time <- unique(data_orig[[time_id]])
  n_ahead <- sum(!(as.character(unique(fframe[[time_id]])) %in% old_time))
  lambda <- model$intercept
  n_labda <- ncol(lambda)
  rho <- model$common_coeff[1]
  n_exp_vars <- length(model$args$exp_var)
  alphas <- model$common_coeff[-1]
  rho_lambda <- c()
  for (i in 0:(n_ahead - 1)) {
    rho_lambda[i + 1] <- rho ^ i
  }
  prediction <- matrix(NA, ncol = n_ahead, nrow = length(lambda))

  # AR(1) model, no exp_vars ---------------------------------------------------
  if (n_exp_vars == 0) {
    T_and_later <- c(
      max(data_orig[[time_id]]),
      head(setdiff(unique(fframe[[time_id]]), old_time), -1)
    )
    for (h in 1:n_ahead) {
      prediction[, h] <- lambda * sum(rho_lambda[1:h]) + rho ^ h * y_T
    }

    # ARX(1) model with horizon T+1 ---------------------------------------------
    # (external variables at T are enough to forecast)
  } else if (n_exp_vars > 0 & n_ahead == 1) {
    T_and_later <- c(
      max(data_orig[[time_id]]),
      head(setdiff(unique(fframe[[time_id]]), old_time), -1)
    )
    Z_impact <- vector("list", n_exp_vars)
    names(Z_impact) <- names(alphas)
    data_last <- dplyr::filter(data_orig, get(time_id) == max(get(time_id)))
    for (a in names(alphas)) {
      Z_impact[[a]] <- alphas[a] * data_last[[a]]
    }
    Z_impact_all_exp_vars <- Reduce("+", Z_impact)
    prediction[, 1] <- lambda + (rho * y_T) + Z_impact_all_exp_vars

    # ARX(1) model with horizon > T+1 -------------------------------------------
    # (new data needed)
  } else if (n_exp_vars > 0 & n_ahead > 1) {
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
    Z_impact <- vector("list", n_exp_vars)
    names(Z_impact) <- names(alphas)
    Z_impact <- lapply(Z_impact, function(x) {
      x <- matrix(NA, ncol = n_ahead, nrow = length(lambda))
    })
    for (a in names(alphas)) {
      for (h in seq(n_ahead)) {
        Z_list <- vector("list", h)
        for (g in 1:h) {
          Z_list[[g]] <- dplyr::filter(newdata2, get(time_id) == T_and_later[g]) %>%
            dplyr::select(a)
        }
        for (z in seq(Z_list)) {
          Z_list[[z]] <- Z_list[[z]] * rev(rho_lambda[1:h])[z]
        }
        sum_h <- Reduce("+", Z_list)
        Z_impact[[a]][, h] <- as.matrix(alphas[a] * sum_h)
      }
    }
    Z_impact_all_exp_vars <- Reduce("+", Z_impact)
    for (h in 1:n_ahead) {
      prediction[, h] <- lambda * sum(rho_lambda[1:h]) + rho ^ h * y_T +
        Z_impact_all_exp_vars[, h]
    }
  }

  # Prepare output ------------------------------------------------------------
  output <- dplyr::select(fframe, unit_id, time_id, model$args$dep_var)
  names(output)[3] <- paste0(names(output)[3], "_TRUE")
  output[[time_id]] <- as.character(output[[time_id]])
  output[[unit_id]] <- as.character(output[[unit_id]])
  fitted_IN <- as.data.frame(model$fitted_values)
  colnames(fitted_IN) <- unique(fframe[[unit_id]])
  rownames(fitted_IN) <- old_time[-1]
  fitted_IN <- pmpp_data(fitted_IN,
    t_dim = "rows",
    var_name = "A"
  )
  fitted_IN$time <- as.character(fitted_IN$time)
  fitted_IN$unit <- as.character(fitted_IN$unit)
  prediction_FC <- as.data.frame(prediction)
  colnames(prediction_FC) <- c(T_and_later[-1], laststamp)
  rownames(prediction_FC) <- unique(fframe[[unit_id]])
  prediction_FC <- pmpp_data(prediction_FC,
    t_dim = "cols",
    var_name = "B"
  )
  prediction_FC$time <- as.character(prediction_FC$time)
  prediction_FC$unit <- as.character(prediction_FC$unit)
  names(prediction_FC)[1:2] <- names(fitted_IN)[1:2] <- names(output)[1:2]
  output <- output %>%
    left_join(fitted_IN, by = c(unit_id, time_id)) %>%
    left_join(prediction_FC, by = c(unit_id, time_id)) %>%
    mutate(temp = ifelse(is.na(A), B, A)) %>%
    dplyr::select(-c(A, B)) %>%
    mutate(window = ifelse(get(time_id) == old_time[1], "EST",
      ifelse(get(time_id) %in% old_time[-1], "IN", "FC")
    ))
  names(output)[4] <- paste0(model$args$dep_var, "_MODEL")
  output[[time_id]] <- factor(output[[time_id]])
  output[[unit_id]] <- factor(output[[unit_id]], levels = 1:model$sample_size["N"])
  output[["window"]] <- factor(output[["window"]])

  return(output)
}
