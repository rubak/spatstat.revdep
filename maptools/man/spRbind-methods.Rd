\name{spRbind-methods}
\docType{methods}
\alias{spRbind-methods}
\alias{spRbind,SpatialPoints,SpatialPoints-method}
\alias{spRbind,SpatialPointsDataFrame,SpatialPointsDataFrame-method}
\alias{spRbind,SpatialLines,SpatialLines-method}
\alias{spRbind,SpatialLinesDataFrame,SpatialLinesDataFrame-method}
\alias{spRbind,SpatialPolygons,SpatialPolygons-method}
\alias{spRbind,SpatialPolygonsDataFrame,SpatialPolygonsDataFrame-method}
\alias{spRbind}
\title{rbind for spatial objects}
\description{
 \code{spRbind} provides rbind-like methods for Spatial*DataFrame objects
}
\section{Methods}{
\describe{

\item{obj = "SpatialPoints", x = "SpatialPoints"}{rbind two SpatialPoints objects}

\item{obj = "SpatialPointsDataFrame", x = "SpatialPointsDataFrame"}{rbind two SpatialPointsDataFrame objects}

\item{obj = "SpatialLines", x = "SpatialLines"}{rbind two SpatialLines objects}

\item{obj = "SpatialLinesDataFrame", x = "SpatialLinesDataFrame"}{rbind two SpatialLinesDataFrame objects}

\item{obj = "SpatialPolygons", x = "SpatialPolygons"}{rbind two SpatialPolygons objects}

\item{obj = "SpatialPolygonsDataFrame", x = "SpatialPolygonsDataFrame"}{rbind two SpatialPolygonsDataFrame objects}
}}

\note{
In addition to the \code{spRbind}-methods, there are also 
\code{rbind}-methods for Spatial* objects. The differences are:

\enumerate{
 \item{\code{spRbind}-methods can bind 2 objects, whereas \code{rbind}-methods 
     can bind multiple object}
 \item{some \code{rbind} can accept objects with duplicated IDs, for all 
     \code{spRbind}-methods these have to be modified explicitly, e.g. by 
     calling \code{\link[sp]{spChFIDs-methods}}}

}}

\author{Roger Bivand}

\seealso{\code{\link[sp]{spChFIDs-methods}}, \code{\link{spCbind-methods}}}

\examples{
xx <- readShapePoly(system.file("shapes/sids.shp", package="maptools")[1], 
  IDvar="FIPSNO", proj4string=CRS("+proj=longlat +ellps=clrk66"))
summary(xx)
xx$FIPSNO
xx1 <- xx[xx$CNTY_ID < 1982,]
xx2 <- xx[xx$CNTY_ID >= 1982,]
xx3 <- spRbind(xx2, xx1)
summary(xx3)
xx3$FIPSNO
}

\keyword{methods}
\keyword{spatial}
