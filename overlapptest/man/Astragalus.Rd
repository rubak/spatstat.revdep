\name{Astragalus}
\alias{Astragalus}
\alias{Sesleria}
\docType{data}
\title{
Astragalus and Sesleria Plants
}
\description{
Oulines of the individuals of \emph{Astragalus sempervirens} and \emph{Sesleria albicans} present in a  2 x 2 m plot in the Pyrenean montains.
}
\usage{
data("Astragalus")
data("Sesleria")
}
\format{
  Each object is a multiple-polygon window with the format \code{\link[spatstat]{owin}} of \pkg{spatstat}. 
}
\examples{
data(Astragalus)
Astragalus
plot(Astragalus)
# total area covered by Astragalus
area.owin(Astragalus)
# number of individual polygons
length(Astragalus$bdry)
# area of each individual 
sapply(Astragalus$bdry, function(x) area.owin(owin(poly=x)))
}
\keyword{datasets}
