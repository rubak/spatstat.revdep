\name{toywin}
\alias{toywin}
\alias{heart}
\alias{bx}
\alias{shp1}
\alias{shp2}
\alias{star}
\docType{data}
\title{
Toy Windows
}
\description{
Synthetic spatial windows for use in testing, simulations and demonstrations.
}
\usage{
data(bx)
data(heart)
data(shp1)
data(shp2)
data(star)
data(toywin)
}
\format{
Each of these is a single closed polygon of class \code{\link[spatstat.geom]{owin}}.
\itemize{
\item{\code{bx} is a box on [-5,5]^2.}
\item{\code{heart} is a heart, professing love for all things \code{\link[spatstat:spatstat-package]{spatstat}}.}
\item{\code{shp1} is shape of mystery.}
\item{\code{shp2} is a slightly more symmetric shape of mystery.}
\item{\code{star} is a star that shines brightly in even non-spatial contexts.}
\item{\code{toywin} is the eponymous toy window used in publications e.g. Davies & Lawson (2019).}
}
}
\details{
These are lazy-loaded so may be called directly by name upon loading of \code{spagmix}.
}
\references{
Davies, T.M. and Lawson, A.B. (2019), An evaluation of likelihood-based bandwidth selectors for spatial and spatiotemporal kernel estimates, \emph{Journal of Statistical Computation and Simulation}, \bold{89} 1131-1152.
}
\author{A.K. Redmond and T.M. Davies}
\examples{
par(mfrow=c(2,3))
plot(bx);axis(1);axis(2)
plot(heart);axis(1);axis(2)
plot(shp1);axis(1);axis(2)
plot(shp2);axis(1);axis(2)
plot(star);axis(1);axis(2)
plot(toywin);axis(1);axis(2)
}
