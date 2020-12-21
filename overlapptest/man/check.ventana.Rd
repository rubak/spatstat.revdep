\name{check.ventana}
\alias{check.ventana}
\title{
Checks for Anticlockwise Vertices
}
\description{
Checks that the vertices of polygons in a multi-polygon \code{\link[spatstat]{owin}} object are listed anticlockwise and, if some are not, tries to correct them.
}
\usage{
check.ventana(ventana)
}
\arguments{
  \item{ventana}{
A multiple-polygon window with the format \code{\link[spatstat]{owin}} of \pkg{spatstat}.
}
}
\details{
This functions should be employed after importing shapefiles into \code{\link[spatstat]{owin}} objetcs (see vignette). Vertices of the individual polygons in the multiple-polygon  \code{\link[spatstat]{owin}} objects should be listed anticlockwise to avoid errors in the computations of area overlapp (clockwise listed polygons represent "holes" in \pkg{spatstat}). This functions checks this and, in case that the vertices of some polygons are listed clockwise, tries to revert their order. 
}
\value{
A multiple-polygon window with the format \code{\link[spatstat]{owin}} of \pkg{spatstat} with the vertices of all polygons listed anticockwise. The order number of the corrected polygons are included in the attribute "malos1".  If there has been any polygon whose vertices have not been corrected, their orden number are included in the attribute "malos2": these polygons should be corrected manually.
}

\author{
Marcelino de la Cruz Rot
}

\examples{
data(Astragalus)
# For illustrative purposes, make the vertices of some individual polygon to be listed clockwise
Astragalus.malo<-Astragalus
Astragalus.malo$bdry[[14]]<-lapply(Astragalus.malo$bdry[[14]], rev)
# check and correct
Astragalus.corrected<-(check.ventana(Astragalus.malo))
attributes(Astragalus.corrected)


}
\keyword{ spatial }
