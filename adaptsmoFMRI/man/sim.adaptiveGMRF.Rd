\name{sim.adaptiveGMRF}
\alias{sim.adaptiveGMRF}
\title{Adaptive GMRF Model for Simulated Data}
\usage{
  sim.adaptiveGMRF(data, hrf, approximate = FALSE, K = 500,
    a = 1, b = 1, c = 1, d = 1, nu = 1, block = 1, burnin =
    1, thin = 1)
}
\arguments{
  \item{data}{simulated fMRI-data, needs to be an array of
  dimension \code{(20 x 20 x T)}.}

  \item{hrf}{haemodynamic response function, needs to be a
  vector of length \code{T}.}

  \item{approximate}{logical, if \code{TRUE} then the
  approximate case is chosen. Default is \code{FALSE}.}

  \item{K}{scalar, length of the MCMC path, hence iteration
  steps.}

  \item{a}{scalar, shape hyperparameter of the
  inverse-gamma distribution of the variance parameter
  (\eqn{\sigma_i^2}).}

  \item{b}{scalar, scale hyperparameter of the inverse
  gamma distribution of the variance parameter
  (\eqn{\sigma_i^2}).}

  \item{c}{scalar, shape hyperparameter of the inverse
  gamma distribution of the precision parameter
  (\eqn{\tau}).}

  \item{d}{scalar, scale hyperparameter of the inverse
  gamma distribution of the precision parameter
  (\eqn{\tau}).}

  \item{nu}{scalar, shape and scale hyperparameter of the
  gamma distribution of the interaction weights
  (\eqn{w_{ij}}).}

  \item{block}{scalar, when \code{approximate==TRUE} then a
  block of weights is updated at a time.}

  \item{burnin}{scalar, defining the first iteration steps
  which should be omitted from MCMC path.}

  \item{thin}{scalar, only every \code{thin} step of MCMC
  path is saved to output.}
}
\description{
  This function estimates the effects of a synthetic
  spatiotemporal data set resembling functional MR Images
  (fMRI), with the method of efficient Markov Chain Monte
  Carlo (MCMC) simulation. The Metropolis Hastings (MH)
  algorithm is used for the non-approximate case and the
  Gibbs sampler for the approximate case.
}
\note{
  This function is solely for one covariate.
}
\examples{
# non-transformed hr-function
T <- 210
seq.length <- T*3
index <- seq(3, T*3, by = 3)
hrf <- rep(c(-0.5, 0.5), each=30, times=ceiling(T/30*1.5))
hrf <- as.matrix(hrf[index])
# get simulated data
data("sim_fmri")
data <- data_simfmri
# execute function
set.seed(111222)
K <- 2
a <- b <- c <- d <- nu <- 1
test.sim.adaptive <- sim.adaptiveGMRF(data, hrf, approximate=TRUE, K,
                                      a, b, c, d, nu)
}
\author{
  Max Hughes
}

