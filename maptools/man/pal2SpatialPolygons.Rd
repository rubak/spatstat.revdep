\name{pal2SpatialPolygons}
\alias{pal2SpatialPolygons}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Making SpatialPolygons objects from RArcInfo input}
\description{
  This function is used in making SpatialPolygons objects from RArcInfo input.
}
\usage{
pal2SpatialPolygons(arc, pal, IDs, dropPoly1=TRUE,
 proj4string=CRS(as.character(NA)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IDs}{Unique character ID values for each output Polygons object; the input IDs can be an integer or character vector with duplicates, where the duplicates will be combined as a single output Polygons object}
  \item{proj4string}{Object of class \code{"CRS"};
	holding a valid proj4 string}
  \item{arc}{Object returned by \code{get.arcdata}}
  \item{pal}{Object returned by \code{get.paldata}}
  \item{dropPoly1}{Should the first polygon in the AVC or e00 data be dropped; the first polygon is typically the compound boundary of the whole dataset, and can be detected by looking at the relative lengths of the list components in the second component of pal, which are the numbers of arcs making up the boundary of each polygon}
}

\value{
  The functions return a SpatialPolygons object
}

\author{ Roger Bivand }

\examples{
nc1 <- readShapePoly(system.file("shapes/sids.shp", package="maptools")[1], ID="FIPS")
plot(nc1)
text(coordinates(nc1), labels=row.names(nc1), cex=0.6)
if(require(maps)){
ncmap <- map("county", "north carolina", fill=TRUE, col="transparent",
 plot=FALSE)
IDs <- sapply(strsplit(ncmap$names, "[,:]"), function(x) x[2])
nc2 <- map2SpatialPolygons(ncmap, IDs)
plot(nc2)
text(coordinates(nc2), labels=row.names(nc2), cex=0.6)
}
#if(require(RArcInfo)) {
#td <- tempdir()
#tmpcover <- paste(td, "nc", sep="/")
#if (!file.exists(tmpcover)) e00toavc(system.file("share/co37_d90.e00",
# package="maptools")[1], tmpcover)
#arc <- get.arcdata(td, "nc")
#pal <- get.paldata(td, "nc")
#pat <- get.tabledata(paste(td, "info", sep="/"), "NC.PAT")
#sapply(pal[[2]], function(x) length(x[[1]]))
#IDs <- paste(pat$ST[-1], pat$CO[-1], sep="")
#nc3 <- pal2SpatialPolygons(arc, pal, IDs=IDs)
#plot(nc3)
#text(coordinates(nc3), labels=row.names(nc3), cex=0.6)
#}
}
\keyword{spatial}
