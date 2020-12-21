\name{rstpoint}
\alias{rstpoint}
\title{
Generate random points in space-time
}
\description{
Generates a random spatiotemporal point pattern containing \eqn{n} independent, identically distributed points with a specified distribution.
}
\usage{
rstpoint(n, f, W = NULL, correction = 1.5, maxpass = 50)
}
\arguments{
  \item{n}{
The number of points to be generated.
}
  \item{f}{
The probability density of the points, an object of class \code{\link{stim}} or \code{\link[sparr:spattemp.density]{stden}}.
}
  \item{W}{
The polygonal \code{\link[spatstat]{owin}} defining the spatial window on which the density is defined. If \code{NULL}, this will be set to the \code{\link[spatstat]{as.polygonal}} version of the pixel images stored in \code{f}. See `Details'.
}
  \item{correction}{
An adjustment to the number of points generated at the initial pass of the internal loop in an effort to minimise the total number of passes required to reach \code{n} points.
}
  \item{maxpass}{
The maximum number of passes allowed before the function exits. If this is reached before \code{n} points are found with respect to the spatiotemporal domain of \code{f}, a warning is issued.
}
}
\details{
This function randomly generates a spatiotemporal point pattern of exactly \code{n} points based on the density function \code{f}. At any given pass, \code{n} * \code{correction} points are generated and rejection sampling is used to accept some of the points; this is repeated until the required number of points is found.

The argument \code{W} is optional, but is useful when the user wants the spatial window of the resulting point pattern to be a corresponding irregular polygon, as opposed to being based on the boundary of a binary image mask (which, when the pixel \code{\link[spatstat]{im}}ages in \code{f} are converted to a polygon directly, gives jagged edges based on the union of the pixels).
}
\value{
An object of class \code{\link[spatstat]{ppp}} containing the \code{n} generated points. The \code{\link[spatstat]{marks}} of the object contain the correspondingly generated observation times.
}
\author{
A.K. Redmond and T.M. Davies
}
\examples{
\donttest{
r1a <- sgmix(cbind(c(0.5,0.5)),vcv=0.01,window=toywin,p0=0.5,p=c(0.5),res=128)
r1b <- sgmix(cbind(c(0.5,0.5),c(0.4,0.6)),vcv=c(0.06,0.015),window=toywin,
             p0=0.1,p=c(0.5,0.4),res=128)
r1c <- sgmix(cbind(c(0.4,0.6)),vcv=c(0.1),window=toywin,p0=0.1,p=c(0.9),res=128)
sts1 <- stkey(start=r1a,
              stop=r1c,
              tlim=c(1,10),
              tres=64,
              window=toywin,
              kf=solist(r1a,r1b),
              kftimes=c(2,6),
              fscale=0.1+0.9*dnorm(seq(-3,3,length=64),mean=0,sd=1))
plot(sts1,sleep=0.1)

Y <- rstpoint(500,sts1,W=toywin,correction=10,maxpass=500)
plot(Y)

require("rgl")
plot3d(Y$x,Y$y,marks(Y))
}
}
