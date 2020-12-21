\name{readShapeLines}
\alias{readShapeLines}
\alias{writeLinesShape}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read arc shape files into SpatialLinesDataFrame objects}
\description{
The use of this function is deprecated and it is not being maintained. Use \code{rgdal::readOGR()} or \code{sf::st_read()} instead - both of these read the coordinate reference system from the input file, while this deprecated function does not. For writing, use \code{rgdal::writeOGR()} or \code{sf::st_write()} instead.

The \code{readShapeLines} function reads data from an arc/line shapefile into a \code{SpatialLinesDataFrame} object; the shapefile may be of type polygon, but for just plotting for example coastlines, a SpatialLines object is sufficient. The \code{writeLinesShape} function writes data from a \code{SpatialLinesDataFrame} object to a shapefile. Note DBF file restrictions in \code{\link[foreign]{write.dbf}}.
}
\usage{
readShapeLines(fn, proj4string=CRS(as.character(NA)), verbose=FALSE,
 repair=FALSE, delete_null_obj=FALSE)
writeLinesShape(x, fn, factor2char = TRUE, max_nchar=254)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fn}{shapefile layer name, when writing omitting the extensions *.shp, *.shx and *.dbf, which are added in the function}
  \item{proj4string}{Object of class \code{CRS}; holding a valid proj4 string}
  \item{verbose}{default FALSE - report type of shapefile and number of shapes}
  \item{repair}{default FALSE: some shapefiles provided by Geolytics Inc. have values of object sizes stored in the *.shx index file that are eight bytes too large, leading the function to try to read past the end of file. If repair=TRUE, an attempt is made to repair the internal values, permitting such files to be read.}
  \item{delete_null_obj}{if TRUE, null geometries will be removed together with their data.frame rows}
  \item{x}{a \code{SpatialLinesDataFrame} object}
  \item{factor2char}{logical, default TRUE, convert factor columns to character}
  \item{max_nchar}{default 254, may be set to a higher limit and passed through to the DBF writer, please see Details in \code{\link[foreign]{write.dbf}}}
}
\details{The shpID values of the shapefile will be used as \code{Lines} ID values; when writing shapefiles, the object data slot row.names are added to the DBF file as column SL\_ID.}

\value{
  a SpatialLinesDataFrame object
}

\author{ Roger Bivand }

\seealso{\code{\link[foreign]{write.dbf}}}


\examples{
xx <- readShapeLines(system.file("shapes/fylk-val.shp", package="maptools")[1],
 proj4string=CRS("+proj=utm +zone=33 +datum=WGS84"))
plot(xx, col="blue")
summary(xx)
xxx <- xx[xx$LENGTH > 30000,]
plot(xxx, col="red", add=TRUE)
tmpfl <- paste(tempdir(), "xxline", sep="/")
writeLinesShape(xxx, tmpfl)
getinfo.shape(paste(tmpfl, ".shp", sep=""))
axx <- readShapeLines(tmpfl, proj4string=CRS("+proj=utm +zone=33 +datum=WGS84"))
plot(xxx, col="black", lwd=4)
plot(axx, col="yellow", lwd=1, add=TRUE)
unlink(paste(tmpfl, ".*", sep=""))
xx <- readShapeLines(system.file("shapes/sids.shp", package="maptools")[1],
 proj4string=CRS("+proj=longlat +datum=NAD27"))
plot(xx, col="blue")
}
\keyword{spatial}

