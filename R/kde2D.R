#' Compute a two-dimensional kernel density estimate
#'
#' @author Michal Oleszak
#'
#' @description The kernel is assumed to be Gaussian. Bandwidth matrix is
#'              diagonal. The two bandwidth parameters are chosen optimally
#'              without ever using/assuming any parametric model for the data
#'              or any "rules of thumb". Unlike many other procedures, this one
#'              is immune to accuracy failures in the estimation of multimodal
#'              densities with widely separated modes. This function in meant to be
#'              the R implementation of the MATLAB \code{kde2d()} function written 
#'              and published by Z. I. Botev at:
#'              \url{http://web.maths.unsw.edu.au/~zdravkobotev/}
#'
#' @references Z. I. Botev, J. F. Grotowski and D. P. Kroese,
#'             "Kernel Density Estimation Via Diffusion",
#'             Annals of Statistics, 2010, Volume 38, Number 5, Pages 2916-2957
#'
#' @param data   N by 2 matrix with the two variables as columns
#' @param n      size of the n by n grid over which the density is computed
#' @param limits limits of the bounding box over which the density is computed;
#'               format: c(lower_Xlim, upper_Xlim, lower_Ylim, upper_Ylim)
#'
#' @return A \code{list} with bandwidth, density and grids for the two dimensions.
#'
#' @importFrom pracma accumarray repmat fzero meshgrid histc numel zeros ones ifft
#' @importFrom stats mvfft
#' @export
kde2D <- function(data, n = 2 ^ 8, limits = NULL) {
  # define auxiliary functions
  ndhist <- function(data, M) {
    bins <- zeros(nrow(data), ncol(data))
    for (i in 1:ncol(data)) {
      bins[, i] <- histc(data[, i], seq(0, 1, 1 / M))$bin
      bins[, i] <- apply(as.matrix(bins[, i]), 1, function(y) min(y, M))
    }
    out <- accumarray(
      bins[all(bins > 0), ], rep((1 / nrow(data)), nrow(data)),
      M * (ones(1, ncol(data)))
    )
    return(out)
  }
  dct2d <- function(data) {
    dct1d <- function(x) {
      x <- rbind(x[seq(1, nrow(x), 2), ], x[seq(nrow(x), 2, -2), ])
      out <- Re(weights * mvfft(x))
      return(out)
    }
    w <- as.matrix(c(1, 2 * (exp(-1i * (1:(nrow(data) - 1)) * pi / (2 * nrow(data))))))
    weights <- w[, ones(1, ncol(data))]
    data <- t(dct1d(t(dct1d(data))))
    return(data)
  }
  idct2d <- function(data) {
    idct1d <- function(x) {
      y <- Re(ifft_mat(weights * x))
      out <- zeros(nrow(data), ncol(data))
      out[seq(1, nrow(data), 2), ] <- y[1:(nrow(data) / 2), ]
      out[seq(2, nrow(data), 2), ] <- y[seq(nrow(data), nrow(data) / 2 + 1, -1), ]
      return(out)
    }
    ifft_mat <- function(x) {
      x <- apply(x, 2, function(y) ifft(y))
      return(x)
    }
    w <- as.matrix(exp(1i * (0:(nrow(data) - 1)) * pi / (2 * nrow(data))))
    weights <- w[, ones(1, ncol(data))]
    out <- idct1d(t(idct1d(data)))
    return(out)
  }
  K <- function(s) {
    if (s == 0) {
      out <- (-1) ^ s / sqrt(2 * pi)
    } else {
      out <- (-1) ^ s * prod(seq(from = 1, to = 2 * s - 1, by = 2)) / sqrt(2 * pi)
    }
    return(out)
  }
  psi <- function(s, time) {
    w <- c((exp(-I * pi ^ 2 * time))) *
      c(cbind(1, 0.5 * ones(1, length(I) - 1)))
    wx <- t(w * (I ^ s[1]))
    wy <- t(w * (I ^ s[2]))
    out <- (-1) ^ sum(s) * (wy %*% A2 %*% t(wx)) * pi ^ (2 * sum(s))
    return(out)
  }
  func <- function(s, t) {
    if (sum(s) <= 4) {
      sum_func <- func(c(s[1] + 1, s[2]), t) + func(c(s[1], s[2] + 1), t)
      const <- (1 + 1 / 2 ^ (sum(s) + 1)) / 3
      time <- c((-2 * const * K(s[1]) * K(s[2]) / N / sum_func) ^ (1 / (2 + sum(s))))
      out <- psi(s, time)
    } else {
      out <- psi(s, t)
    }
    return(out)
  }
  evolve <- function(t) {
    sum_func <- func(c(0, 2), t) + func(c(2, 0), t) + 2 * func(c(1, 1), t)
    time <- (2 * pi * N * sum_func) ^ (-1 / 3)
    out <- (t - time) / time
    return(c(out, time))
  }
  subtract_evolve <- function(t) {
    return(t - evolve(t))
  }

  # compute the kernel density estimate
  N <- nrow(data)
  # round up n to the next power of 2
  n <- 2 ^ ceiling(log2(n))
  # define default limits
  if (is.null(limits)) {
    max <- c(max(data[, 1]), max(data[, 2]))
    min <- c(min(data[, 1]), min(data[, 2]))
    range <- max - min
    max_xy <- max + range / 4
    min_xy <- min - range / 4
  } else {
    max_xy <- c(limits[2], limits[4])
    min_xy <- c(limits[1], limits[3])
  }
  # scale data
  scaling <- max_xy - min_xy
  data_trans <- (data - repmat(min_xy, N, 1)) / repmat(scaling, N, 1)
  # bin the data uniformly using regular grid
  data_init <- ndhist(data_trans, n)
  # discrete cosine transform of initial data
  data_DCT <- dct2d(data_init)
  # compute the optimal bandwidth ^ 2
  I <- (0:(n - 1)) ^ 2
  A2 <- data_DCT ^ 2
  t_star <- fzero(subtract_evolve, c(0, 0.1))[[1]]
  p_02 <- func(c(0, 2), t_star)
  p_20 <- func(c(2, 0), t_star)
  p_11 <- func(c(1, 1), t_star)
  t_y <- (p_02 ^ (3 / 4) / (4 * pi * N * p_20 ^ (3 / 4) * (p_11 + sqrt(p_20 * p_02)))) ^ (1 / 3)
  t_x <- (p_20 ^ (3 / 4) / (4 * pi * N * p_02 ^ (3 / 4) * (p_11 + sqrt(p_20 * p_02)))) ^ (1 / 3)
  # smooth the discrete cosine transform of initial data using t_star
  data_DCT_smoothed <- (t(t(exp(-(0:(n - 1)) ^ 2 * pi ^ 2 * c(t_x) / 2))) %*%
    t(exp(-(0:(n - 1)) ^ 2 * pi ^ 2 * c(t_y) / 2))) * data_DCT
  # apply the inverse discrete cosine transform
  density <- idct2d(data_DCT_smoothed) *
    (numel(data_DCT_smoothed) / prod(scaling))
  grid <- meshgrid(
    seq(min_xy[1], max_xy[1], by = scaling[1] / (n - 1)),
    seq(min_xy[2], max_xy[2], by = scaling[2] / (n - 1))
  )
  # compute the bandwidth
  bw <- sqrt(cbind(t_x, t_y)) * scaling

  return(list(
    "bandwidth" = bw, "density" = density,
    "X" = grid[[1]], "Y" = grid[[2]]
  ))
}
