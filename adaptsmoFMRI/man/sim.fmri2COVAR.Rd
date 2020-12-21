\name{sim.fmri2COVAR}
\alias{sim.fmri2COVAR}
\title{Simulate FMRI Data}
\usage{
  sim.fmri2COVAR(hrf, beta.Var1, beta.Var2)
}
\arguments{
  \item{hrf}{haemodynamic response function, needs to be a
  vector of length \code{T}.}

  \item{beta.Var1}{scalar, defines the height of the
  activated area, in form of a cylinder of the first grid.}

  \item{beta.Var2}{scalar, defines the height of the
  activated area, in form of a cylinder of the second
  grid.}
}
\description{
  This function returns the synthetic spatiotemporal data
  set resembling functional MR Images (fMRI) data.
}
\details{
  The returned data is simulated on a 20 x 20 grid.
}
\note{
  This function is solely for two covariates.
}
\examples{
# non-transformed hr-function
T <- 180
seq.length <- T*3
index <- seq(3, T*3, by = 3)
vis <- rep(c(-0.5, 0.5), each=30, times=ceiling(T/30*1.5))
vis <- as.matrix(vis[index])
aud <- rep(c(-0.5, 0.5), each=45, times=ceiling(T/30*1.5))
aud <- as.matrix(aud[index])
hrf <- cbind(vis,aud)
# define height of activation area
beta.Var1 <- beta.Var2 <- 3
# use function to obtain fmri data
data <- sim.fmri2COVAR(hrf, beta.Var1, beta.Var2)$fmri
}
\author{
  Max Hughes
}

