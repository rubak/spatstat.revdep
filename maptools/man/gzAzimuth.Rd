\name{gzAzimuth}
\alias{gzAzimuth}
\alias{trackAzimuth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find azimuth for geographical coordinates}
\description{
  The function finds azimuth values for geographical coordinates given as decimal degrees from the \code{from} coordinates to the \code{to} coordinate. In function \code{trackAzimuth}, the azimuth values are found between successive rows of the input coordinate matrix.
}
\usage{
gzAzimuth(from, to, type = "snyder_sphere")
trackAzimuth(track, type="snyder_sphere")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{from}{a two column matrix of geographical coordinates given as decimal degrees (longitude first)}
  \item{track}{a two column matrix of geographical coordinates given as decimal degrees (longitude first)}
  \item{to}{a one row, two column matrix or two element vector of geographical coordinates given as decimal degrees (longitude first)}
  \item{type}{default is \code{"snyder_sphere"}, otherwise \code{"abdali"}; the results should be identical with slightly less trigonometry in \code{"abdali"}}
}
\details{
  The azimuth is calculated on the sphere, using the formulae given by Snyder (1987, p. 30) and Abdali (1997, p. 17). The examples use data taken from Abdali (p. 17--18). There is a very interesting discussion of the centrality of azimuth-finding in the development of mathematics and mathematical geography in Abdali's paper. Among others, al-Khwarizmi was an important contributor. As Abdali puts it, "This is a veritable who's who of medieval science" (p. 3).
}
\value{
  values in decimal degrees - zero is North - of the azimuth from the \code{from} coordinates to the \code{to} coordinate.
}
\references{Snyder JP (1987) Map projections - a working manual, USGS Professional Paper 1395; Abdali SK (1997) "The Correct Qibla", formerly at http://patriot.net/users/abdali/ftp/qibla.pdf}
\author{Roger Bivand, with contributions by Sebastian Luque}

\examples{
name <- c("Mecca", "Anchorage", "Washington")
long <- c(39.823333, -149.883333, -77.0166667)
lat <- c(21.423333, 61.2166667, 38.9)
x <- cbind(long, lat)
row.names(x) <- name
crib <- c(-9.098363, 56.575960)
r1 <- gzAzimuth(x[2:3,], x[1,])
r1
all.equal(r1, crib)
r2 <- gzAzimuth(x[2:3,], x[1,], type="abdali")
r2
all.equal(r2, crib)
trackAzimuth(x)
}
\keyword{spatial}
