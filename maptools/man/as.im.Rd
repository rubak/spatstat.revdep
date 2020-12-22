\name{as.im}
\alias{as.im.SpatialGridDataFrame}
\alias{coerce,SpatialGridDataFrame,im-method}
\alias{as.SpatialGridDataFrame.im}
\alias{coerce,im,SpatialGridDataFrame-method}
\alias{as.im.RasterLayer}

\title{Coercion between sp objects and spatstat im objects}
\description{
Functions to convert between \pkg{spatstat}s raster format \code{im} and \pkg{sp}s \code{SpatialGridDataFrame} as well as one-way conversion from \pkg{raster}s \code{rasterLayer} to \code{im}. S4-style \code{as()} coercion can be used between \code{im} and \code{SpatialGridDataFrame} objects.
}
\section{Methods}{
  \describe{
\item{coerce}{\code{signature(from = "SpatialGridDataFrame", to = "im")}}
\item{coerce}{\code{signature(from = "im", to = "SpatialGridDataFrame")}}
}
}

\usage{
as.im.SpatialGridDataFrame(from)
as.SpatialGridDataFrame.im(from)
as.im.RasterLayer(from, factor.col.name = NULL)
}
\arguments{
  \item{from}{object to coerce from}
  \item{factor.col.name}{column name of \code{levels(from)} to be treated as a factor; if NULL defaults to last column of \code{from}. Ignored if \code{from} is not a raster with factor values.}
}
\details{A \code{SpatialGridDataFrame} object may contain several columns of 
data such that several values are associated with each grid cell. In contrast an
\code{im} object can only contain a single variable value for each cell. In
\code{as.im.SpatialGridDataFrame()} the first data column is used. To convert 
another column to \code{im} format simply extract this column first as shown in
the example below.}

\section{Warning}{
  In \pkg{spatstat} all spatial objects are assumed to be planar. This means
  that \pkg{spatstat} is not designed to work directly with geographic
  (longitude and latitude) coordinates. If a \pkg{sp} object is declared to
  have geographic (unprojected) coordinates \pkg{maptools} refuses to convert
  directly to \pkg{spatstat} format. Rather, these should be projected first
  using e.g. \code{\link[sp]{spTransform}}. If you know what you are doing, and
  really want to force coercion, you can overwrite the \code{proj4string} of the
  \pkg{sp} object with \code{NA}, \code{proj4string(x) <- CRS(NA)}, which
  will fool the system to think that the data is in local planar coordinates.
  This is probably not a good idea!
}

\seealso{
  Other converters between \pkg{sp} and \pkg{spatstat}:
  \code{\link{as.ppp.SpatialPoints}},
  \code{\link{as.psp.SpatialLines}},
  \code{\link{as.owin.SpatialPolygons}},
  \code{\link{as.SpatialPolygons.tess}}.
}

\author{Edzer Pebesma \email{edzer.pebesma@uni-muenster.de}, Roger Bivand}
\examples{
run <- FALSE
if (require(spatstat.linnet, quietly=TRUE)) run <- TRUE
if (run) {
## Extract an example SpatialGridDataFrame and plot it
data(meuse.grid) # A data.frame
gridded(meuse.grid) = ~x+y # Now a SpatialPixelsDataFrame
fullgrid(meuse.grid) <- TRUE # Finally a SpatialGridDataFrame
mg_dist <- meuse.grid["dist"] # A SpatialGridDataFrame with a single column
image(mg_dist, axes=TRUE)
}
if (run) {
## Convert to im format and plot it
mg_im <- as(mg_dist, "im")
plot(mg_im)
}
if (run) {
## Convert back to SpatialGridDataFrame and plot it again
mg2 <- as.SpatialGridDataFrame.im(mg_im)
image(mg2, axes=TRUE)
}
run <- run && require(raster, quietly=TRUE)
if (run) {
## Convert SpatialGridDataFrame -> RasterLayer -> im and plot it
r <- as(mg2, "RasterLayer")
r_im <- as.im.RasterLayer(r)
plot(r_im)
}
if (run) {
rr <- raster(nrow=2, ncol=3)
values(rr) <- 1:6
rr <- as.factor(rr)
rrr <- rr
f <- levels(rrr)[[1]]
f$name <- c("vector", "of", "values")
f$name2 <- letters[1:6]
levels(rrr) <- f
print(levels(rrr))
}
if (run) {
iii <- as.im.RasterLayer(rrr)
plot(iii)
}
if (run) {
iv <- as.im.RasterLayer(rrr, factor.col.name = "name")
plot(iv)
}
if (run) {
}
}
\keyword{spatial}
