\name{symbolsInPolys}
\alias{symbolsInPolys}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Place grids of points over polygons}
\description{
  Place grids of points over polygons with chosen density and/or symbols (suggested by Michael Wolf).
}
\usage{
symbolsInPolys(pl, dens, symb = "+", compatible = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pl}{an object of class SpatialPolygons or SpatialPolygonsDataFrame}
  \item{dens}{number of symbol plotting points per unit area; either a single numerical value for all polygons, or a numeric vector the same length as pl with values for each polygon}
  \item{symb}{plotting symbol; either a single value for all polygons, or a vector the same length as pl with values for each polygon}
  \item{compatible}{what to return, if TRUE a a list of matrices of point coordinates, one matrix for each member of pl, with a \code{symb} attribute, if false a SpatialPointsDataFrame with a \code{symb} column}
}
\details{
  The dots are placed in a grid pattern with the number of points per polygon being polygon area times density (number of dots not guaranteed to be the same as the count). When the polygon is made up of more than one part, the dots will be placed in proportion to the relative areas of the clockwise rings (anticlockwise are taken as holes). From maptools release 0.5-2, correction is made for holes in the placing of the dots, but depends on hole values being correctly set, which they often are not. 
}

\value{
The function returns a list of matrices of point coordinates, one matrix for each member of pl; each matrix has a \code{symb} attribute that can be used for setting the \code{pch} argument for plotting. If the count of points for the given density and polygon area is zero, the list element is NULL, and can be tested when plotting - see the examples.
}

\note{Extension to plot pixmaps at the plotting points using addlogo() from the pixmap package is left as an exercise for the user.}

\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link[sp]{spsample}}}

\examples{
nc_SP <- readShapePoly(system.file("shapes/sids.shp", package="maptools")[1],
 proj4string=CRS("+proj=longlat +ellps=clrk66"))
\dontrun{
pls <- slot(nc_SP, "polygons")
pls_new <- lapply(pls, checkPolygonsHoles)
nc_SP <- SpatialPolygonsDataFrame(SpatialPolygons(pls_new,
 proj4string=CRS(proj4string(nc_SP))), data=as(nc_SP, "data.frame"))
}
symbs <- c("-", "+", "x")
np <- sapply(slot(nc_SP, "polygons"), function(x) length(slot(x, "Polygons")))
try1 <- symbolsInPolys(nc_SP, 100, symb=symbs[np])
plot(nc_SP, axes=TRUE)
plot(try1, add=TRUE, pch=as.character(try1$symb))
}
\keyword{spatial}

