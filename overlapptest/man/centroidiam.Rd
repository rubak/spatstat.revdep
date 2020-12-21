\name{centroidiam}
\alias{centroidiam}
\title{
Compute Centroids and Diameters
}
\description{
Computes the centroid and diameter of each individual polygon in a multi-polygon \code{\link[spatstat]{owin}} object.
}
\usage{
centroidiam(ventana1)
}

\arguments{
  \item{ventana1}{
A multiple-polygon window with the format \code{\link[spatstat]{owin}} of \pkg{spatstat}
}
}
\details{
Iteratively applies the functions \code{\link[spatstat]{centroid.owin}} and \code{\link[spatstat]{diameter}} of \pkg{spatstat} to each polygon in the multipolygon \code{\link[spatstat]{owin}} and computes its centroid and its diameter.
}
\value{
 \item{diams }{Vector of diameters}
  \item{centroids }{Matrix with the coordinates of the centroids}

}

\author{
Marcelino de la Cruz Rot
}

\examples{
data(Astragalus)
X<-centroidiam(Astragalus)
X$centroids
X$diams
}

\keyword{ spatial }
