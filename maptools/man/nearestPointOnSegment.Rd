\name{nearestPointOnSegment}
\alias{nearestPointOnSegment}
\title{
Get the nearest point on a segment to a given point
}
\description{
This function calculates the coordinates of and the distance to the nearest 
point on a segment to a given point. This function does not work with 
geographic coordinates.
}
\usage{
nearestPointOnSegment(s, p)
}
\arguments{
  \item{s}{
A matrix representing the coordinates of the segment. The matrix has 2x2 
dimension where each row represents one of the end points. 
}
  \item{p}{
A vector representing the X and Y coordinates of the point.
}
}
\value{
A vector with three numeric values representing X and Y coordinates of the 
nearest point on a segment to a given point as well as the distance between
both points.
}
\references{
The function was ported to R based on this code: http://pastebin.com/n9rUuGRh
}
\author{
German Carrillo 
}
\seealso{
\code{\link{nearestPointOnLine}}, \code{\link{snapPointsToLines}} 
}
\examples{
segment = cbind(c(1,2),c(1,1.5))
point = c(1.2,1.5)
nearestPointOnSegment(segment, point)
}
\keyword{ spatial }
