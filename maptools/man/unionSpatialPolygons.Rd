\name{unionSpatialPolygons}
\alias{unionSpatialPolygons}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Aggregate Polygons in a SpatialPolygons object}
\description{
  The function aggregates Polygons in a SpatialPolygons object, according to the IDs vector specifying which input Polygons belong to which output Polygons; internal boundaries are dissolved using the rgeos package \code{gUnaryUnion} function. If the rgeos package is not available, and if the gpclib package is available and the user confirms that its restrictive license conditions are met, its \code{union} function will be used.
}
\usage{
unionSpatialPolygons(SpP, IDs, threshold=NULL, avoidGEOS=FALSE, avoidUnaryUnion=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SpP}{A SpatialPolygons object as defined in package sp}
  \item{IDs}{A vector defining the output Polygons objects, equal in length to the length of the polygons slot of SpRs; it may be character, integer, or factor (try table(factor(IDs)) for a sanity check). It may contain NA values for input objects not included in the union}
  \item{threshold}{if not NULL, an area measure below which slivers will be discarded (some polygons have non-identical boundaries, for instance along rivers, generating slivers on union which are artefacts, not real sub-polygons)}
  \item{avoidGEOS}{default FALSE; if TRUE force use of gpclib even when GEOS is available}
  \item{avoidUnaryUnion}{avoid gUnaryUnion if it is available; not relevant before GEOS 3.3.0}
}

\value{
  Returns an aggregated SpatialPolygons object named with the aggregated IDs values in their sorting order; see the ID values of the output object to view the order.
}

\section{Warning}{When using GEOS Unary Union, it has been found that some polygons are not dissolved correctly when the absolute values of the coordinates are very small. No work-around is available at present.}

\author{Roger Bivand}
%\seealso{\code{\link[rgeos]{gUnaryUnion}}}
\examples{
if (rgeosStatus()) {
nc1 <- readShapePoly(system.file("shapes/sids.shp", package="maptools")[1],
 proj4string=CRS("+proj=longlat +datum=NAD27"))
lps <- coordinates(nc1)
ID <- cut(lps[,1], quantile(lps[,1]), include.lowest=TRUE)
reg4 <- unionSpatialPolygons(nc1, ID)
row.names(reg4)
}
}
\keyword{spatial}
