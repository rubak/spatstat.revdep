\name{sim.weightsplot}
\alias{sim.weightsplot}
\title{Plot Function for Weights (Simulated Data)}
\usage{
  sim.weightsplot(weights, nei, nx, ny, thresh = 0.05, ...)
}
\arguments{
  \item{weights}{matrix, containing MCMC-output the of
  posterior estimates of weights.}

  \item{nei}{matrix, locations of weights in precision
  matrix.}

  \item{nx}{scalar, number of pixels in x-direction.}

  \item{ny}{scalar, number of pixels in y-direction.}

  \item{thresh}{scalar, defining the threshold to which the
  median of the weights smaller than this threshold should
  be plotted.}

  \item{\dots}{graphical parameters for \code{image} can
  also be passed on as arguments to this function.}
}
\description{
  This function plots the weights interacting between
  estimated effects for each pixel.
}
\note{
  This function is solely for MCMC-outputs on simulated
  data.
}
\author{
  Max Hughes
}

