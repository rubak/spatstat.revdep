\name{spCbind-methods}
\docType{methods}
\alias{spCbind-methods}
\alias{spCbind,SpatialPointsDataFrame,data.frame-method}
\alias{spCbind,SpatialPointsDataFrame,vector-method}
\alias{spCbind,SpatialLinesDataFrame,data.frame-method}
\alias{spCbind,SpatialLinesDataFrame,vector-method}
\alias{spCbind,SpatialPolygonsDataFrame,data.frame-method}
\alias{spCbind,SpatialPolygonsDataFrame,vector-method}
\title{cbind for spatial objects}
\alias{spCbind}
\description{
 \code{spCbind} provides cbind-like methods for Spatial*DataFrame objects in addition to the \code{$}, \code{[<-} and \code{[[<-} methods already available.
}
\section{Methods}{
\describe{

\item{obj = "SpatialPointsDataFrame", x = "data.frame"}{cbind a data frame to the data slot of a SpatialPointsDataFrame object}

\item{obj = "SpatialPointsDataFrame", x = "vector"}{cbind a vector to the data slot of a SpatialPointsDataFrame object}

\item{obj = "SpatialLinesDataFrame", x = "data.frame"}{cbind a data frame to the data slot of a SpatialLinesDataFrame object; the data frame argument must have row names set to the Lines ID values, and should be re-ordered first by matching against a shared key column}

\item{obj = "SpatialLinesDataFrame", x = "vector"}{cbind a vector to the data slot of a SpatialLinesDataFrame object}

\item{obj = "SpatialPolygonsDataFrame", x = "data.frame"}{cbind a data frame to the data slot of a SpatialPolygonsDataFrame object; the data frame argument must have row names set to the Polygons ID values, and should be re-ordered first by matching against a shared key column}

\item{obj = "SpatialPolygonsDataFrame", x = "vector"}{cbind a vector to the data slot of a SpatialPolygonsDataFrame object}
}}

\author{Roger Bivand}

\seealso{\code{\link[sp]{spChFIDs-methods}}, \code{\link{spRbind-methods}}}

\examples{
xx <- readShapePoly(system.file("shapes/sids.shp", package="maptools")[1], 
  IDvar="FIPSNO", proj4string=CRS("+proj=longlat +ellps=clrk66"))
library(foreign)
xtra <- read.dbf(system.file("share/nc_xtra.dbf", package="maptools")[1])
o <- match(xx$CNTY_ID, xtra$CNTY_ID)
xtra1 <- xtra[o,]
row.names(xtra1) <- xx$FIPSNO
xx1 <- spCbind(xx, xtra1)
names(xx1)
identical(xx1$CNTY_ID, xx1$CNTY_ID.1)
}

\keyword{methods}
\keyword{spatial}
