\name{rgmix}
\alias{rgmix}
\title{Random bivariate Gaussian mixture density generation}
\description{
Generates a pixel image of a bivariate normal mixture density observed on a bounded window using a specified number of contributing densities with randomly selected means and variance-covariance matrices.
}
\usage{
rgmix(N, window, v = 4, S = NULL, extras = FALSE, ...)
}
\arguments{
\item{N}{
  The number of Gaussian components to generate for the mixture.
}

\item{window}{
  An object of class \code{\link[spatstat.geom]{owin}} giving the observational window on which the mixture density is defined.
}

\item{v}{
  The degrees of freedom for the inverse-Wishart distribution of the variance-covariance matrices (must be at least 4). The default value of 4 ensures the generated covariance matrices are centered on \code{S}.
}
\item{S}{
  A symmetric, positive-definite \eqn{2 \times 2} scale matrix for the inverse-Wishart distribution of the variance-covariance matrices.
}
\item{extras}{
  A logical value indicating whether, in addition to returning the pixel \code{\link[spatstat.geom]{im}}age of the final mixture density, to also return the randomly realised mean locations and corresponding variance-covariance matrices. See `Value'.
}
\item{...}{
  Additional arguments to be passed to \code{\link{sgmix}}. See `Details'.
}
}
\details{
This function creates and returns a bivariate Gaussian mixture density on a bounded \code{window} based on \code{N} randomly generated mean locations and corresponding randomly generated variance-covariance matrices. First, the \code{N} mean locations are generated based on a uniform distribution over the spatial \code{window}. Each location is then associated with a covariance matrix generated from an inverse-Wishart distribution with \code{v} degrees of freedom and scale matrix \code{S}.

Once the above steps are completed, the function calls \code{\link{sgmix}} with the chosen mean and covariance matrices, thereby creating the Gaussian mixture. Resolution and other aspects of this call can be controlled by using \code{...}, passing the contents internally to \code{\link{sgmix}}. By default, all generated Gaussian components have equal weight in contributing to the final mixture density. The user can alter this by passing \code{p0} and \code{p} to the \code{...}, though should take care that the length of \code{p} is \code{N}, and that \code{p0} and \code{p} sum to 1, as outlined in the documentation for \code{\link{sgmix}}.
}

\value{
If \code{extras = FALSE} (default), then a pixel \code{\link[spatstat.geom]{im}}age of the final mixture density. If \code{extras = TRUE}, a list is returned with members \code{f} (the pixel \code{\link[spatstat.geom]{im}}age of the final mixture density); \code{mn} (a \eqn{2 \times} \code{N} matrix with each column giving the mean location of each of the \code{N} Gaussian bumps); and \code{vcv} (a \eqn{2 \times 2 \times} \code{N} array with layers giving the covariance matrices associated with the means in the columns of \code{mn}).
}

\author{A.K. Redmond and T.M. Davies}

\examples{
set.seed(321)
dens1 <- rgmix(7,window=toywin)
plot(dens1)

set.seed(456)
dens2 <- rgmix(7,window=toywin)
plot(dens2)

# Explicitly return details of generated means and covariances
set.seed(321)
dens1.detailed <- rgmix(7,window=toywin,extras=TRUE)
dens1.detailed$f
dens1.detailed$mn
dens1.detailed$vcv

# Set underlying uniform proportion and compare with dens2 from above
set.seed(456)
dens2.wunif <- rgmix(7,window=toywin,p0=0.3)
par(mfrow=c(1,2))
plot(rpoint(500,dens2))
plot(rpoint(500,dens2.wunif))

# Explicitly setting scale matrix for inverse-wishart generation of covariances
dens3 <- rgmix(3,window=toywin,S=matrix(c(0.025,-0.004,-0.004,0.02),2))
plot(dens3)
}
