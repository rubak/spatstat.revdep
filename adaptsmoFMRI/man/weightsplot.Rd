\name{weightsplot}
\alias{weightsplot}
\title{Plot Function for Weights (Real Data)}
\usage{
  weightsplot(weights, nei, nx, ny, coord, thresh = 0.1,
    ...)
}
\arguments{
  \item{weights}{matrix, containing MCMC-output the of
  posterior estimates of weights.}

  \item{nei,}{matrix, locations of weights in precision
  matrix.}

  \item{nx,}{scalar, number of pixels in x-direction.}

  \item{ny,}{scalar, number of pixels in y-direction.}

  \item{coord,}{matrix, coordinates of pixels.}

  \item{thresh,}{scalar, defining the threshold to which
  the median of the weights smaller than this threshold
  should be plotted.}

  \item{\dots,}{graphical parameters for \code{image} can
  also be passed on as arguments to this function.}
}
\description{
  This function plots the weights interacting between
  estimated effects for each pixel.
}
\note{
  This function is solely for MCMC-outputs on real data.
}
\author{
  Max Hughes
}

