\name{stgmix}
\alias{stgmix}
\alias{stim}
\title{Trivariate Gaussian mixture density generation}
\description{Generates a pixel image array of a specified trivariate normal mixture density observed on a bounded window in space and time.}
\usage{
stgmix(mean, vcv, window, tlim, p0 = 0, p = NULL, sres = 128, tres = sres, int = 1)
}
\arguments{
  \item{mean}{
A \eqn{3 x N} matrix specifying the means of each of \eqn{N} contributing normal densities; each component in the order of (x-coord, y-coord, time-coord).
}
  \item{vcv}{
A \eqn{3 x 3 x N} array specifying the variance-covariance matrices of each contributing density.
}
  \item{window}{
An object of class \code{\link[spatstat]{owin}} giving the spatial observational window on which the mixture density is defined.
}
  \item{tlim}{
A vector of length 2 giving the boundaries of the time interval on which the mixture density is defined.
}
  \item{p0}{
The proportion of uniform density that contributes to the final mixture (default is 0).
}
  \item{p}{
  A numeric vector of the \eqn{N} proportions for each contributing density (default is equal proportions for each density, after subtracting \code{p0}). Together, \code{p0} and \code{p} must sum to exactly 1.
}
  \item{sres}{
The spatial resolution (number of pixels) along each side of the spatial grid (default is 128).
}
  \item{tres}{
The temporal resolution (default is to equate with \code{sres}).
}
  \item{int}{
  A positive numeric value for post-hoc rescaling of the density (useful if the user wishes to return a spatiotemporal \emph{intensity} function). Defaults to 1 for no change in scaling.
}
}


\details{
This function creates a 3D array of a density function made up of a mixture of \eqn{N} trivariate normals with the interpretation of a continuous probability density function in space-time. As such, each component is restricted to conserve mass over a 3D region specified by a fixed polygonal \code{window} in space, stretched over defined temporal limits (\code{tlim}). A warning will appear if less than 1\% of the integral of each Gaussian bump is inside this observational spatiotemporal polyhedron.
}

\value{
An object of class \code{stim} giving the trivariate density. This is a list with six components:
\item{a}{The \code{sres} \eqn{x} \code{sres} \eqn{x} \code{tres} array of the specified density.}
\item{v}{A pixel \code{\link[spatstat]{im}}age version of \code{a}, provided as a \code{\link[spatstat]{solist}} of length \code{tres}, with each member being the spatial \code{\link[spatstat]{im}}age slice of the 3D density at each of the time-coordinate values.}
\item{xcol}{Grid coordinates in the spatial x-axis (corresponds to each spatial \code{\link[spatstat]{im}}age in \code{v}).}
\item{yrow}{Grid coordinates in the spatial y-axis (corresponds to each spatial \code{\link[spatstat]{im}}age in \code{v}).}
\item{tlay}{Grid coordinates in the temporal axis (corresponds to the order of the spatial \code{\link[spatstat]{im}}ages in \code{v}).}
\item{W}{A copy of \code{window}, the spatial \code{\link[spatstat]{owin}} upon which the density is defined.}
}

\author{A.K. Redmond and T.M. Davies}

\examples{
\donttest{
require("abind")
m1 <- c(0.3,0.3,2)
m2 <- c(0.5,0.8,8)
m3 <- c(0.7,0.6,7)
v1 <- diag(c(0.01^2,0.01^2,1))
v2 <- diag(c(0.005,0.005,0.5))
v3 <- diag(c(0.005,0.005,0.5))
stg1 <- stgmix(mean=cbind(m1,m2,m3),
               vcv=abind(v1,v2,v3,along=3),
               window=toywin,tlim=c(1,10),
               p0=0.1,tres=64)
plot(stg1,log=TRUE)


mn <- matrix(c(0,0,0,-2,1,4,1,-2,8),nrow=3)
vr <- array(c(1,0,0,0,1,0,0,0,1,1,0,0.5,0,1,0,0.5,0,3,1,0,0,0,2,0,0,0,1),
            dim=c(3,3,3))
stg2 <- stgmix(mean=mn,vcv=vr,window=shp1,
               tlim=c(0,10),tres=50)
plot(stg2,fix.range=TRUE,sleep=0.1)
}
}
