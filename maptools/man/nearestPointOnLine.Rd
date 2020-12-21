\name{nearestPointOnLine}
\alias{nearestPointOnLine}
\title{
Get the nearest point on a line to a given point
}
\description{
This function calculates the coordinates of the nearest point on a line to a 
given point. This function does not work with geographic coordinates.
}
\usage{
nearestPointOnLine(coordsLine, coordsPoint)
}
\arguments{
  \item{coordsLine}{
Matrix with coordinates of line vertices. Each row represents a vertex. 
}
  \item{coordsPoint}{
A vector representing the X and Y coordinates of the point.
}
}
\value{
Vector with the X and Y coordinates of the nearest point on a line to the 
given point.
}
\author{
German Carrillo
}
\seealso{
\code{\link{nearestPointOnSegment}}, \code{\link{snapPointsToLines}} 
}
\examples{
coordsLine = cbind(c(1,2,3),c(3,2,2))
coordsPoint = c(1.2,1.5)
nearestPointOnLine(coordsLine, coordsPoint)
}
\keyword{ spatial }
