#' Transform a single variable in the matrix format into the long panel format
#'
#' @author Michal Oleszak
#'
#' @description This function transforms a matrix of data with cross-sectional
#'              and time dimensions in rows and columns or columns and rows into
#'              a panel-structured, 3-column data frame
#'
#' @param indata   matrix with  a single variable
#' @param t_dim    character string, one of: 'cols', 'rows';
#'                 whether time dimension in indata is across columns or rows
#' @param var_name character string; name of the variable in indata
#'
#' @return A \code{data.frame} with 3 columns: unit, time and variable's values.
#' @export
#'
#' @examples
#' set.seed(1)
#' matrix_var <- matrix(rnorm(100), nrow = 20)
#' panel_var <- pmpp_data(matrix_var)
pmpp_data <- function(indata, t_dim = "cols", var_name = "Y") {
  if (!is.matrix(indata)) {
    indata <- as.matrix(indata)
  }
  if (t_dim == "cols") {
    N <- nrow(indata)
    T <- ncol(indata)
    Y <- c()
    for (i in 1:N) {
      Y <- rbind(Y, as.matrix(indata[i, ], ncol = 1))
    }
    if (is.null(colnames(indata))) {
      colnames(indata) <- 1:T
    }
    if (is.null(rownames(indata))) {
      rownames(indata) <- 1:N
    }
    out <- data.frame(
      rep(rownames(indata), each = T),
      rep(colnames(indata), N), as.numeric(Y)
    )
  } else if (t_dim == "rows") {
    N <- ncol(indata)
    T <- nrow(indata)
    Y <- c()
    for (i in 1:N) {
      Y <- rbind(Y, as.matrix(indata[, i], ncol = 1))
    }
    if (is.null(colnames(indata))) {
      colnames(indata) <- 1:N
    }
    if (is.null(rownames(indata))) {
      rownames(indata) <- 1:T
    }
    out <- data.frame(
      rep(colnames(indata), each = T),
      rep(rownames(indata), N), as.numeric(Y)
    )
  }
  colnames(out) <- c("unit", "time", var_name)
  return(out)
}
