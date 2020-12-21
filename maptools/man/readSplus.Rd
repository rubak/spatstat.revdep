\name{readSplus}
\alias{readSplus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read exported WinBUGS maps}
\description{
  The function permits an exported WinBUGS map to be read into an \pkg{sp} package class \code{SpatialPolygons} object.
}
\usage{
readSplus(file, proj4string = CRS(as.character(NA)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{name of file}
  \item{proj4string}{Object of class '"CRS"'; holding a valid proj4 string}
}

\value{
  \code{readSplus} returns a SpatialPolygons object
}
\references{\url{http://www.mrc-bsu.cam.ac.uk/wp-content/uploads/geobugs12manual.pdf}}
\author{Virgilio Gomez Rubio <Virgilio.Gomez@uclm.es>}
\note{In the example, taken from the GeoBUGS manual, the smaller part of area1 has a counter-clockwise ring direction in the data, while other rings are clockwise. This implies that it is a hole, and does not get filled. Errant holes may be filled using \code{\link{checkPolygonsHoles}}. The region labels are stored in the \code{ID} slots of the \code{Polygons} objects.}
\seealso{\code{\link{map2SpatialPolygons}}}
\examples{
if (rgeosStatus()) {
geobugs <- readSplus(system.file("share/Splus.map", package="maptools"))
plot(geobugs, axes=TRUE, col=1:3)
row.names(geobugs)
pls <- slot(geobugs, "polygons")
sapply(pls, function(i) sapply(slot(i, "Polygons"), slot, "hole"))
pls1 <- lapply(pls, checkPolygonsHoles)
sapply(pls1, function(i) sapply(slot(i, "Polygons"), slot, "hole"))
plot(SpatialPolygons(pls1), axes=TRUE, col=1:3)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

