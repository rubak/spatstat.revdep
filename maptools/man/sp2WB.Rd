\encoding{latin1}
\name{sp2WB}
\alias{sp2WB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Export SpatialPolygons object as S-Plus map for WinBUGS}
\description{
  The function exports an sp SpatialPolygons object into a S-Plus map
format to be import by WinBUGS.
}
\usage{
sp2WB(map, filename, Xscale = 1, Yscale = Xscale, plotorder = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{map}{a SpatialPolygons object}
  \item{filename}{file where output is written}
  \item{Xscale, Yscale}{scales to be written in the output file}
  \item{plotorder}{default=FALSE, if TRUE, export polygons in plotting order}
}

\references{\url{http://www.mrc-bsu.cam.ac.uk/wp-content/uploads/geobugs12manual.pdf}}
\author{Virgilio Gómez Rubio, partly derived from earlier code by Thomas Jagger}


\examples{
xx <- readShapePoly(system.file("shapes/sids.shp", package="maptools")[1],
 IDvar="FIPSNO", proj4string=CRS("+proj=longlat +ellps=clrk66"))
plot(xx, border="blue", axes=TRUE, las=1)
tf <- tempfile()
sp2WB(as(xx, "SpatialPolygons"), filename=tf)
xxx <- readSplus(tf, proj4string=CRS("+proj=longlat +ellps=clrk66"))
all.equal(xxx, as(xx, "SpatialPolygons"), tolerance=.Machine$double.eps^(1/4),
 check.attributes=FALSE)
\dontrun{
x <- readAsciiGrid(system.file("grids/test.ag", package="maptools")[1])
xp <- as(x, "SpatialPixelsDataFrame")
pp <- as(xp, "SpatialPolygons")
td <- tempdir()
sp2WB(pp, filename=file.path(td, "test.map"))
}
}
\keyword{spatial}
