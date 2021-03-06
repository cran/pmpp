% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kde2D.R
\name{kde2D}
\alias{kde2D}
\title{Compute a two-dimensional kernel density estimate}
\usage{
kde2D(data, n = 2^8, limits = NULL)
}
\arguments{
\item{data}{N by 2 matrix with the two variables as columns}

\item{n}{size of the n by n grid over which the density is computed}

\item{limits}{limits of the bounding box over which the density is computed;
format: c(lower_Xlim, upper_Xlim, lower_Ylim, upper_Ylim)}
}
\value{
A \code{list} with bandwidth, density and grids for the two dimensions.
}
\description{
The kernel is assumed to be Gaussian. Bandwidth matrix is
             diagonal. The two bandwidth parameters are chosen optimally
             without ever using/assuming any parametric model for the data
             or any "rules of thumb". Unlike many other procedures, this one
             is immune to accuracy failures in the estimation of multimodal
             densities with widely separated modes. This function in meant to be
             the R implementation of the MATLAB \code{kde2d()} function written 
             and published by Z. I. Botev at:
             \url{http://web.maths.unsw.edu.au/~zdravkobotev/}
}
\references{
Z. I. Botev, J. F. Grotowski and D. P. Kroese,
            "Kernel Density Estimation Via Diffusion",
            Annals of Statistics, 2010, Volume 38, Number 5, Pages 2916-2957
}
\author{
Michal Oleszak
}
