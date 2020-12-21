\name{SpatialLinesMidPoints}
\alias{SpatialLinesMidPoints}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Line midpoints}
\description{
The function onverts SpatialLinesDataFrame to SpatialPointsDataFrame with points at the midpoints of the line segments.
}
\usage{
SpatialLinesMidPoints(sldf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sldf}{A SpatialLines or SpatialLinesDataFrame object}
}
\details{
The function builds a SpatialPointsDataFrame from the midpoints of Line objects belonging to Lines objects in an object inheriting from a Spatial Lines object. The output data slot contains an index variable showing which Lines object the midpoints belong to.
}

\value{
  A SpatialPointsDataFrame object created from the input object.
}

\author{Jonathan Callahan, modified by Roger Bivand}

\examples{
xx <- readShapeLines(system.file("shapes/fylk-val.shp", package="maptools")[1],
 proj4string=CRS("+proj=utm +zone=33 +datum=WGS84"))
plot(xx, col="blue")
spdf <- SpatialLinesMidPoints(xx)
plot(spdf, col="orange", add=TRUE)
}

\keyword{spatial}

