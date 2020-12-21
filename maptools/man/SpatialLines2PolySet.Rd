\name{SpatialLines2PolySet}
\alias{SpatialLines2PolySet}
\alias{SpatialPolygons2PolySet}
\alias{PolySet2SpatialPolygons}
\alias{PolySet2SpatialLines}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert sp line and polygon objects to PBSmapping PolySet objects}
\description{
  Functions \code{SpatialLines2PolySet} and \code{SpatialPolygons2PolySet} convert objects of sp classes to \code{PolySet} class objects as defined in the PBSmapping package, and \code{PolySet2SpatialLines} and \code{PolySet2SpatialPolygons} convert in the opposite direction.
}
\usage{
SpatialLines2PolySet(SL)
SpatialPolygons2PolySet(SpP)
PolySet2SpatialLines(PS)
PolySet2SpatialPolygons(PS, close_polys=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SL}{a SpatialLines object as defined in the sp package}
  \item{SpP}{a SpatialPolygons object as defined in the sp package}
  \item{PS}{a PolySet object}
  \item{close_polys}{should polygons be closed if open}
}

\value{
  PolySet objects as defined in the PBSmapping package
}

\author{Roger Bivand and Andrew Niccolai}


\seealso{\code{\link[PBSmapping]{PolySet}}, \code{\link{MapGen2SL}}}

\examples{
if(require(PBSmapping) && require(maps)) {
nor_coast_lines <- map("world", interior=FALSE, plot=FALSE, xlim=c(4,32),
 ylim=c(58,72))
nor_coast_lines <- pruneMap(nor_coast_lines, xlim=c(4,32), ylim=c(58,72))
nor_coast_lines_sp <- map2SpatialLines(nor_coast_lines,
 proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
nor_coast_lines_PS <- SpatialLines2PolySet(nor_coast_lines_sp)
summary(nor_coast_lines_PS)
plotLines(nor_coast_lines_PS)
o3 <- PolySet2SpatialLines(nor_coast_lines_PS)
plot(o3, axes=TRUE)
nor_coast_poly <- map("world", "norway", fill=TRUE, col="transparent",
 plot=FALSE, ylim=c(58,72))
IDs <- sapply(strsplit(nor_coast_poly$names, ":"), function(x) x[1])
nor_coast_poly_sp <- map2SpatialPolygons(nor_coast_poly, IDs=IDs,
 proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
nor_coast_poly_PS <- SpatialPolygons2PolySet(nor_coast_poly_sp)
summary(nor_coast_poly_PS)
plotPolys(nor_coast_poly_PS)
o1 <- PolySet2SpatialPolygons(nor_coast_poly_PS)
plot(o1, axes=TRUE)
}
}
\keyword{spatial}

