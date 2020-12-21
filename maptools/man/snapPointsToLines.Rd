\name{snapPointsToLines}
\alias{snapPointsToLines}
\title{
Snap a set of points to a set of lines
}
\description{
This function snaps a set of points to a set of lines based on the minimum 
distance of each point to any of the lines. This function does not work with 
geographic coordinates.
}
\usage{
snapPointsToLines(points, lines, maxDist=NA, withAttrs = TRUE, idField=NA)
}
\arguments{
  \item{points}{
An object of the class SpatialPoints or SpatialPointsDataFrame.
}
  \item{lines}{
An object of the class SpatialLines or SpatialLinesDataFrame.
}
  \item{maxDist}{
Numeric value for establishing a maximum distance to avoid snapping points that 
are farther apart; its default value is NA.
}
  \item{withAttrs}{
Boolean value for preserving (TRUE) or getting rid (FALSE) of the original point 
attributes. Default: TRUE. This parameter is optional. 
}
  \item{idField}{
A string specifying the field which contains each line's id. This id will 
be transferred to the snapped points data set to distinguish the line which each
point was snapped to.
}
}
\value{
SpatialPointsDataFrame object as defined by the R package 'sp'. This object 
contains the snapped points, therefore all of them lie on the lines.
}
\author{
German Carrillo and Ethan Plunkett
}
\seealso{
\code{\link{nearestPointOnSegment}}, \code{\link{nearestPointOnLine}}, 
\code{\link{sp}}
}
\examples{
# From the sp vignette
l1 = cbind(c(1,2,3),c(3,2,2))
l1a = cbind(l1[,1]+.05,l1[,2]+.05)
l2 = cbind(c(1,2,3),c(1,1.5,1))
Sl1 = Line(l1)
Sl1a = Line(l1a)
Sl2 = Line(l2)
S1 = Lines(list(Sl1, Sl1a), ID="a")
S2 = Lines(list(Sl2), ID="b")
Sl = SpatialLines(list(S1,S2))
df = data.frame(z = c(1,2), row.names=sapply(slot(Sl, "lines"), function(x) slot(x, "ID")))
Sldf = SpatialLinesDataFrame(Sl, data = df)

xc = c(1.2,1.5,2.5)
yc = c(1.5,2.2,1.6)
Spoints = SpatialPoints(cbind(xc, yc))

if (rgeosStatus()) snapPointsToLines(Spoints, Sldf, maxDist=0.4)
}
\keyword{ spatial }
