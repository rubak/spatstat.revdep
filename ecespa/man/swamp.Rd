\name{swamp}
\alias{swamp}
\docType{data}
\title{ Tree Species in a Swamp Forest }
\description{
  Locations and botanical classification of trees in a plot in the Savannah River. Locations are given in metres,
  rounded to the nearest 0.1 metre.  The data come from  a 1-ha (200 m x 50 m) plot in the Savannah River Site, South Carolina, USA. 
The 734 mapped stems included 156 Carolina ash (\emph{Fraxinus caroliniana}), 
215 Water tupelo (\emph{Nyssa aquatica}), 205 Swamp tupelo (\emph{Nyssa sylvatica}), 
98 Bald cypress (\emph{Taxodium distichum}) and 60 stems of 8 additional species.
Although the plots were set up by Bill Good and their spatial patterns described in Good
and Whipple(1982), the plots have been maintained and resampled by Rebecca Sharitz and her colleagues of the Savannah River Ecology Laboratory.
There are slightly different versions of the Good plot data. Every time the plots are
resampled, some errors are corrected. This is mostly a concern for the
biologists. The different versions are very similar; they are all very good
examples of a marked spatial point pattern.

}
\usage{data(swamp)}
\format{
  A data frame with 734 observations on the following 3 variables.
  \describe{
    \item{\code{x}}{Cartesian x-coordinate of tree}
    \item{\code{y}}{Cartesian y-coordinate of tree}
    \item{\code{sp}}{a factor with levels indicating the species
 of each tree:}
 \tabular{ll}{
     \code{FX}\tab Carolina ash (\emph{Fraxinus caroliniana})\cr
	 \code{NS}\tab Swamp tupelo (\emph{Nyssa sylvatica})\cr
	 \code{NX}\tab Water tupelo (\emph{Nyssa aquatica})\cr
	 \code{TD}\tab Bald cypress (\emph{Taxodium distichum})\cr
	 \code{OT}\tab Other species\cr
	}
  }
}

\references{
Dixon, P.M. 2002. Nearest-neighbor contingency table analysis of spatial segregation for several species.
  \emph{Ecoscience}, 9 (2): 142-151. \url{https://doi.org/10.1080/11956860.2002.11682700}. 

Good, , B. J. & Whipple, S.A. 1982. Tree spatial patterns: South Carolina bottomland
and swamp forest. \emph{Bulletin of the Torrey Botanical Club}, 109: 529-536.

Jones et al. 1994. Tree population dynamics in seven South Carolina mixed-species
forests. \emph{Bulletin of the Torrey Botanical Club}, 121:360-368.

}
\examples{
data(swamp)
plot(swamp$x,swamp$y, col=as.numeric(swamp$sp),pch=19,
		 xlab="",ylab="",main="Swamp forest")
}
\keyword{datasets}
