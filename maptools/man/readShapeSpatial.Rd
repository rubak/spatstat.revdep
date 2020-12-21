\name{readShapeSpatial}
\alias{readShapeSpatial}
\alias{writeSpatialShape}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read shape files into Spatial*DataFrame objects}
\description{
The use of this function is deprecated and it is not being maintained. Use \code{rgdal::readOGR()} or \code{sf::st_read()} instead - both of these read the coordinate reference system from the input file, while this deprecated function does not. For writing, use \code{rgdal::writeOGR()} or \code{sf::st_write()} instead.

The \code{readShapeSpatial} reads data from a shapefile into a \code{Spatial*DataFrame} object. The \code{writeSpatialShape} function writes data from a \code{Spatial*DataFrame} object to a shapefile. Note DBF file restrictions in \code{\link[foreign]{write.dbf}}.
}
\usage{
readShapeSpatial(fn, proj4string=CRS(as.character(NA)), 
	verbose=FALSE, repair=FALSE, IDvar=NULL, force_ring=FALSE, 
	delete_null_obj=FALSE, retrieve_ABS_null=FALSE)
writeSpatialShape(x, fn, factor2char = TRUE, max_nchar=254)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fn}{shapefile layer name, when writing omitting the extensions *.shp, *.shx and *.dbf, which are added in the function}
  \item{proj4string}{Object of class \code{CRS}; holding a valid proj4 string}
  \item{verbose}{default FALSE - report type of shapefile and number of shapes}
  \item{repair}{default FALSE: some shapefiles provided by Geolytics Inc. have values of object sizes stored in the *.shx index file that are eight bytes too large, leading the function to try to read past the end of file. If repair=TRUE, an attempt is made to repair the internal values, permitting such files to be read.}
  \item{IDvar}{a character string: the name of a column in the shapefile DBF containing the ID values of the shapes - the values will be converted to a character vector (Polygons only)}
  \item{force_ring}{if TRUE, close unclosed input rings (Polygons only)}
  \item{delete_null_obj}{if TRUE, null geometries inserted by ABS will be removed together with their data.frame rows (Polygons and Lines)}
  \item{retrieve_ABS_null}{default FALSE, if TRUE and delete\_null\_obj also 
TRUE, the function will return a data frame containing the data from any 
null geometries inserted by ABS (Polygons only)}
  \item{x}{a vector data \code{Spatial*DataFrame} object}
  \item{factor2char}{logical, default TRUE, convert factor columns to character}
  \item{max_nchar}{default 254, may be set to a higher limit and passed through to the DBF writer, please see Details in \code{\link[foreign]{write.dbf}}}
}
\details{If no IDvar argument is given, the shpID values of the shapefile will be used as \code{Polygons} ID values; when writing shapefiles, the object data slot row.names are added to the DBF file as column SP\_ID.}
\value{
   a Spatial*DataFrame object of a class corresponding to the input shapefile
}

\author{ Roger Bivand }

\seealso{\code{\link[foreign]{write.dbf}}}

\examples{
library(maptools)
xx <- readShapeSpatial(system.file("shapes/sids.shp", package="maptools")[1],
 IDvar="FIPSNO", proj4string=CRS("+proj=longlat +ellps=clrk66"))
summary(xx)
xxx <- xx[xx$SID74 < 2,]
tmpfl <- paste(tempdir(), "xxpoly", sep="/")
writeSpatialShape(xxx, tmpfl)
getinfo.shape(paste(tmpfl, ".shp", sep=""))
unlink(paste(tmpfl, ".*", sep=""))
xx <- readShapeSpatial(system.file("shapes/fylk-val.shp",
 package="maptools")[1], proj4string=CRS("+proj=utm +zone=33 +datum=WGS84"))
summary(xx)
xxx <- xx[xx$LENGTH > 30000,]
plot(xxx, col="red", add=TRUE)
tmpfl <- paste(tempdir(), "xxline", sep="/")
writeSpatialShape(xxx, tmpfl)
getinfo.shape(paste(tmpfl, ".shp", sep=""))
unlink(paste(tmpfl, ".*", sep=""))
xx <- readShapeSpatial(system.file("shapes/baltim.shp", package="maptools")[1])
summary(xx)
xxx <- xx[xx$PRICE < 40,]
tmpfl <- paste(tempdir(), "xxpts", sep="/")
writeSpatialShape(xxx, tmpfl)
getinfo.shape(paste(tmpfl, ".shp", sep=""))
unlink(paste(tmpfl, ".*", sep=""))
}
\keyword{spatial}
