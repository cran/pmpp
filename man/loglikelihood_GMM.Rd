% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pmpp_aux_funcs.R
\name{loglikelihood_GMM}
\alias{loglikelihood_GMM}
\title{Produce negative log-likelihood in the GMM case}
\usage{
loglikelihood_GMM(theta, rho_GMMpar, alpha_GMMpar, sigma2_GMMpar, n_alpha,
  X_mat, Y_mat, Z_mat, W, T, N, aux_Y0)
}
\arguments{
\item{theta}{vector of homogeneous parameters}

\item{rho_GMMpar}{lagged dependent variables coefficient estimates from the GMM}

\item{alpha_GMMpar}{external variables coefficient estimates from the GMM}

\item{sigma2_GMMpar}{variance of the shocks estimated using GMM residuals}

\item{n_alpha}{number of external variables}

\item{X_mat}{lagged dependent variable matrix}

\item{Y_mat}{dependent variable matrix}

\item{Z_mat}{external variable matrix}

\item{W}{cross-sectionally invariant variables - not used now}

\item{T}{time dimension of the data}

\item{N}{cross-sectional dimension of the data}

\item{aux_Y0}{auxiliary matrix with initial observations of the dependent variable}
}
\description{
Produce negative log-likelihood in the GMM case
}
