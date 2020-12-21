\name{sgmix}
\alias{sgmix}
\title{Bivariate Gaussian mixture density generation}
\description{Generates a pixel image of a specified bivariate normal mixture density observed on a bounded window.}
\usage{
sgmix(mean, vcv, window, p0 = 0, p = NULL, resolution = 128, int = 1)
}
\arguments{
  \item{mean}{
  A \eqn{2 \times N} matrix specifying the means of each of \eqn{N} contributing normal densities.
  }
  \item{vcv}{
  Either a \eqn{2 \times 2 \times N} array specifying the variance-covariance matrices of each contributing density, or a numeric vector of length \eqn{N} giving the isotropic \bold{standard deviations} of each contributing density. An error is thrown if the function encounters anything but a symmetric, positive-definite covariation specification for each component.
  }
  \item{window}{
  An object of class \code{\link[spatstat]{owin}} giving the observational window on which the mixture density is defined.
  }
  \item{p0}{
  The proportion of uniform density that contributes to the mixture (default is 0).
  }
  \item{p}{
  A numeric vector of the \eqn{N} proportions for each contributing density (default is equal proportions for each density, after subtracting \code{p0}). Together, \code{p0} and \code{p} must sum to exactly 1.
  }
  \item{resolution}{
  The number of pixels along each side of the grid for the pixel image (default is 128).
  }
  \item{int}{
  A positive numeric value for post-hoc rescaling of the density (useful if the user wishes to return an \emph{intensity} function). Defaults to 1 for no change in scaling.
  }
}

\details{
This function generates a pixel \code{\link[spatstat]{im}}age of a 2D density function made of a mixture of \eqn{N} bivariate normals; each component is restricted to conserve probability mass over a bounded subset of the plane. A warning will appear if less than 1\% of the integral of each Gaussian bump is inside the observational window.}

\value{
An object of class \code{\link[spatstat]{im}} giving the mixture density.
}

\author{
A.K. Redmond
}

\examples{
# Example using isotropic standard deviations
m1 <- c(0.4,0.5)
m2 <- c(0.2,0.7)
s1 <- 0.1
s2 <- 0.025
dens1 <- sgmix(mean=cbind(m1,m2),vcv=c(s1,s2),window=toywin,p0=0.3,p=c(0.5,0.2))

plot(dens1,log=TRUE)
pts1 <- rpoint(200,dens1) # generate random points via spatstat::rpoint
points(pts1)


# Example using full covariance matrices
mn <- cbind(c(0.25,0.8),c(0.31,0.82),c(0.43,0.64),c(0.63,0.62),c(0.21,0.26))
v1 <- matrix(c(0.0023,-0.0009,-0.0009,0.002),2)
v2 <- matrix(c(0.0016,0.0015,0.0015,0.004),2)
v3 <- matrix(c(0.0007,0.0004,0.0004,0.0011),2)
v4 <- matrix(c(0.0023,-0.0041,-0.0041,0.0099),2)
v5 <- matrix(c(0.0013,0.0011,0.0011,0.0014),2)
vr <- array(NA,dim=c(2,2,5))
for(i in 1:5) vr[,,i] <- get(paste("v",i,sep=""))
dens2 <- sgmix(mean=mn,vcv=vr,window=toywin,p0=0.1)

plot(dens2,log=TRUE)
pts2 <- rpoint(200,dens2)
points(pts2)
}
