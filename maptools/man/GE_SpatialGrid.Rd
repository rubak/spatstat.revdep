\name{GE_SpatialGrid}
\alias{GE_SpatialGrid}
\alias{Sobj_SpatialGrid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create SpatialGrid for PNG output to GE}
\description{
  The function sets up metadata in the form of a SpatialGrid object for defining the size and placing of a PNG image overlay in Google Earth. The internal function \code{Sobj\_SpatialGrid} can also be called to build a grid for arbitrary Spatial* objects.
}
\usage{
GE_SpatialGrid(obj, asp = NA, maxPixels = 600)
Sobj_SpatialGrid(obj, asp=1, maxDim=100, n=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{a Spatial* object}
  \item{asp}{if NA, will be set to the latitude corrected value}
  \item{maxPixels}{the maximum dimension of the output PNG}
  \item{maxDim}{the maximum dimension of the output grid; ignored if \code{n} not NULL}
  \item{n}{if not NULL, the minimum number of cells in the returned grid}
}
\details{
  The function is used together with \code{kmlOverlay} to wrap around the opening of a PNG graphics device, plotting code, and the closing of the device. The computed values take account of the adjustment of the actual data bounding box to an integer number of rows and columns in the image file. 

The approach may be used as an alternative to writing PNG files from SpatialGrid and SpatialPixel objects in \pkg{rgdal} using \code{writeGDAL}, and to writing KML files using \code{writeOGR} for vector data objects. The output PNG files are likely to be very much smaller than large vector data KML files, and hinder the retrieval of exact positional information.

Note that the geometries should be in geographical coordinates with datum WGS84 for export to KML.
}
\value{
  returns an S3 object of class \code{GE_SG} with components:
  \item{height}{Integer raster height for png call}
  \item{width}{Integer raster width for png call}
  \item{SG}{a SpatialGrid object with the grid topology of the output PNG}
  \item{asp}{the aspect value used}
  \item{xlim}{xlim taken from SG}
  \item{ylim}{ylim taken from SG}
}
\author{Duncan Golicher, David Forrest and Roger Bivand}
\seealso{\code{\link{kmlOverlay}}}
\examples{
opt_exask <- options(example.ask=FALSE)
qk <- SpatialPointsDataFrame(quakes[, c(2:1)], quakes)
summary(Sobj_SpatialGrid(qk)$SG)
t2 <- Sobj_SpatialGrid(qk, n=10000)$SG
summary(t2)
prod(slot(slot(t2, "grid"), "cells.dim"))
proj4string(qk) <- CRS("+proj=longlat +ellps=WGS84")
tf <- tempfile()
SGqk <- GE_SpatialGrid(qk)
png(file=paste(tf, ".png", sep=""), width=SGqk$width, height=SGqk$height,
  bg="transparent")
par(mar=c(0,0,0,0), xaxs="i", yaxs="i")
plot(qk, xlim=SGqk$xlim, ylim=SGqk$ylim, setParUsrBB=TRUE)
dev.off()
kmlOverlay(SGqk, paste(tf, ".kml", sep=""), paste(tf, ".png", sep=""))
\dontrun{
qk0 <- quakes
qk0$long <- ifelse(qk0$long <= 180, qk0$long, qk0$long-360)
qk0a <- SpatialPointsDataFrame(qk0[, c(2:1)], qk0)
proj4string(qk0a) <- CRS("+proj=longlat +ellps=WGS84")
# writeOGR(qk0a, paste(tf, "v.kml", sep=""), "Quakes", "KML")
# system(paste("googleearth ", tf, ".kml", sep=""))
}
options(example.ask=opt_exask)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

