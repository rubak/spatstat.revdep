\name{sim.fmri}
\alias{sim.fmri}
\title{Simulate FMRI Data}
\usage{
  sim.fmri(hrf, beta)
}
\arguments{
  \item{hrf}{haemodynamic response function, needs to be a
  vector of length \code{T}.}

  \item{beta}{scalar, defines the height of the activated
  area, in form of a cylinder.}
}
\description{
  This function returns the synthetic spatiotemporal data
  set resembling functional MR Images (fMRI) data.
}
\details{
  The returned data is simulated on a 20 x 20 grid.
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
# define height of activation area
beta <- 3
# use function to obtain fmri data
data <- sim.fmri(hrf, beta)$fmri
}
\author{
  Max Hughes
}

