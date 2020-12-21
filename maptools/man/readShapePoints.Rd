\name{readShapePoints}
\alias{readShapePoints}
\alias{writePointsShape}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read points shape files into SpatialPointsDataFrame objects}
\description{
The use of this function is deprecated and it is not being maintained. Use \code{rgdal::readOGR()} or \code{sf::st_read()} instead - both of these read the coordinate reference system from the input file, while this deprecated function does not.For writing, use \code{rgdal::writeOGR()} or \code{sf::st_write()} instead.

The \code{readShapePoints} reads data from a points shapefile into a \code{SpatialPointsDataFrame} object. The \code{writePointsShape} function writes data from a \code{SpatialPointsDataFrame} object to a shapefile. Both reading and writing can be carried out for 2D and 3D point coordinates. Note DBF file restrictions in \code{\link[foreign]{write.dbf}}.
}
\usage{
readShapePoints(fn, proj4string = CRS(as.character(NA)), verbose = FALSE,
 repair=FALSE)
writePointsShape(x, fn, factor2char = TRUE, max_nchar=254)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fn}{shapefile layer name, when writing omitting the extensions *.shp, *.shx and *.dbf, which are added in the function}
  \item{proj4string}{Object of class \code{CRS}; holding a valid proj4 string}
  \item{verbose}{default FALSE - report type of shapefile and number of shapes}
  \item{repair}{default FALSE: some shapefiles provided by Geolytics Inc. have values of object sizes stored in the *.shx index file that are eight bytes too large, leading the function to try to read past the end of file. If repair=TRUE, an attempt is made to repair the internal values, permitting such files to be read.}
  \item{x}{a \code{SpatialPointsDataFrame} object}
  \item{factor2char}{logical, default TRUE, convert factor columns to character}
  \item{max_nchar}{default 254, may be set to a higher limit and passed through to the DBF writer, please see Details in \code{\link[foreign]{write.dbf}}}
}

\value{
  a SpatialPointsDataFrame object
}

\author{ Roger Bivand }

\seealso{\code{\link[foreign]{write.dbf}}}

\examples{
library(maptools)
xx <- readShapePoints(system.file("shapes/baltim.shp", package="maptools")[1])
plot(xx)
summary(xx)
xxx <- xx[xx$PRICE < 40,]
tmpfl <- paste(tempdir(), "xxpts", sep="/")
writePointsShape(xxx, tmpfl)
getinfo.shape(paste(tmpfl, ".shp", sep=""))
axx <- readShapePoints(tmpfl)
plot(axx, col="red", add=TRUE)
unlink(paste(tmpfl, ".*", sep=""))
xx <- readShapePoints(system.file("shapes/pointZ.shp", package="maptools")[1])
dimensions(xx)
plot(xx)
summary(xx)
}
\keyword{spatial}

