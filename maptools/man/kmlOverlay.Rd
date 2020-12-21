\name{kmlOverlay}
\alias{kmlOverlay}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create and write KML file for PNG image overlay}
\description{
 The function is used to create and write a KML file for a PNG image overlay for Google Earth.
}
\usage{
kmlOverlay(obj, kmlfile = NULL, imagefile = NULL, name = "R image")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{a \code{GE_SG} object from \code{GE_SpatialGrid}}
  \item{kmlfile}{if not NULL the name of the kml file to be written}
  \item{imagefile}{the name of the PNG file containing the image - this should be either relative (same directory as kml file) or abosolute (fully qualified)}
  \item{name}{the name used to describe the image overlay in GE}
}
\details{
    The function is used together with \code{GE_SpatialGrid} to wrap around the opening of a PNG graphics device, plotting code, and the closing of the device. The computed values take account of the adjustment of the actual data bounding box to an integer number of rows and columns in the image file. 

The approach may be used as an alternative to writing PNG files from SpatialGrid and SpatialPixel objects in \pkg{rgdal} using \code{writeGDAL}, and to writing KML files using \code{writeOGR} for vector data objects. The output PNG files are likely to be very much smaller than large vector data KML files, and hinder the retrieval of exact positional information.

Note that the geometries should be in geographical coordinates with datum WGS84.
}
\value{
x is a character vector containing the generated lines of the kml file
}

\author{Duncan Golicher, David Forrest and Roger Bivand}

\seealso{\code{\link{GE_SpatialGrid}}}
\examples{
opt_exask <- options(example.ask=FALSE)
qk <- SpatialPointsDataFrame(quakes[, c(2:1)], quakes)
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
#library(rgdal)
#qk0 <- quakes
#qk0$long <- ifelse(qk0$long <= 180, qk0$long, qk0$long-360)
#qk0a <- SpatialPointsDataFrame(qk0[, c(2:1)], qk0)
#proj4string(qk0a) <- CRS("+proj=longlat +ellps=WGS84")
#writeOGR(qk0a, paste(tf, "v.kml", sep=""), "Quakes", "KML")
#system(paste("googleearth ", tf, ".kml", sep=""))
}
options(example.ask=opt_exask)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

