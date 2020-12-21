\name{rotawin}
\alias{rotawin}
\title{
Rotate Individual Polygons
}
\description{
Randomly rotates individual polygons around their centroids.
}
\usage{
rotawin(ventana)
}

\arguments{
  \item{ventana}{
A multiple-polygon window with the format \code{\link[spatstat.geom]{owin}} of \pkg{spatstat}.
}
}
\details{
\code{rotawin} applies an independent random rotation to each of the polygons in a multiple-polygon \code{\link[spatstat.geom]{owin}} object.
}
\value{
\code{rotawin} returns the original \code{\link[spatstat.geom]{owin}} object with the individual polygons randomly rotated.
}
\author{
Marcelino de la Cruz Rot
}

\examples{
data(Astragalus)
plot(Astragalus)
plot(rotawin(Astragalus), add=TRUE, border=2)

}
\keyword{ spatial }

