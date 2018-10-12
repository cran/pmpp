#' Summary method for objects of class \code{pmpp}.
#'
#' @param object object of class \code{pmpp}, as returned by \code{pmpp()}
#' @param file a connection, or a character string naming the file to print to
#' @param ... other parameters passed further
#'
#' @return A summary object for class \code{pmpp}.
#' @export
#' @importFrom moments skewness kurtosis
#' @importFrom stats shapiro.test IQR median
summary.pmpp <- function(object, file = "", ...) {
  if (!inherits(object, "pmpp")) {
    stop("Non-convenient object, PMPP object required.")
  }
  cat(paste("\nCall:\n"), file = file)
  cat(paste(deparse(object$call), "\n"), file = file, append = TRUE)
  cat(paste("\nElapsed time:\n"), file = file)
  cat(paste((paste0(object$elapsed_time, "\n"))))
  cat(paste("\nCommon coefficients:\n"), file = file)
  print(object$common_coeff)
  cat(paste("\nIntercept:\n"), file = file)
  int <- data.frame(
    c("Mean:", "Median:", "IQR:", "Skewness:", "Kurtosis:"),
    round(c(
      mean(object$intercept), median(object$intercept),
      IQR(object$intercept), skewness(c(object$intercept)),
      kurtosis(c(object$intercept), na.rm = TRUE)
    ), 3)
  )
  names(int) <- NULL
  print(int)
  cat(paste("\nResiduals:\n"), file = file)
  err <- data.frame(
    c("RMSE (in-sample):", "MAPE (in-sample):", "Residuals normality p-value (H0: normal):"),
    round(c(object$RMSE, object$MAPE, shapiro.test(object$residuals)$p.value), 3)
  )
  names(err) <- NULL
  print(err)
}
