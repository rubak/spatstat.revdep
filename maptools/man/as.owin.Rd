\name{as.owin}
\alias{as.owin.SpatialGridDataFrame}
\alias{coerce,SpatialGridDataFrame,owin-method}
\alias{as.owin.SpatialPixelsDataFrame}
\alias{coerce,SpatialPixelsDataFrame,owin-method}
\alias{as.owin.SpatialPolygons}
\alias{coerce,SpatialPolygons,owin-method}
\alias{as.SpatialPolygons.owin}
\alias{coerce,owin,SpatialPolygons-method}

\title{Coercion between sp objects and spatstat owin objects}
\description{
Functions to convert between \pkg{spatstat}s observation window (\code{owin})
format and various \pkg{sp} formats. S4-style \code{as()} coercion can be used
as well.
}
\section{Methods}{
  \describe{
\item{coerce}{\code{signature(from = "SpatialPolygons", to = "owin")}}
\item{coerce}{\code{signature(from = "SpatialPixelsDataFrame", to = "owin")}}
\item{coerce}{\code{signature(from = "SpatialGridDataFrame", to = "owin")}}
\item{coerce}{\code{signature(from = "owin", to = "SpatialPolygons")}}
}
}

\usage{
as.owin.SpatialPolygons(W, \dots, fatal)
as.owin.SpatialGridDataFrame(W, \dots, fatal)
as.owin.SpatialPixelsDataFrame(W, \dots, fatal)
as.SpatialPolygons.owin(x)
}
\arguments{
  \item{W}{\code{SpatialPolygons} object to coerce to \code{owin}}
  \item{x}{\code{owin} object to coerce to \code{SpatialPolygons} format}
  \item{\dots}{ignored}
  \item{fatal}{formal coercion argument; ignored}
}

\section{Warning}{
  In \pkg{spatstat} all spatial objects are assumed to be planar. This means
  that \pkg{spatstat} is not designed to work directly with geographic
  (longitude and latitude) coordinates. If a \pkg{sp} object is declared to
  have geographic (unprojected) coordinates \pkg{maptools} refuses to convert
  directly to \pkg{spatstat} format. Rather, these should be projected first
  using e.g. \code{\link[sp]{spTransform}}. If you know what you are doing, and
  really want to force coercion, you can overwrite the \code{proj4string} of the
  \pkg{sp} object with an empty string, \code{proj4string(x) <- ""}, which
  will fool the system to think that the data is in local planar coordinates.
  This is probably not a good idea!
}


\details{An observation window in \pkg{spatstat} defines a planar region. It is
typically used to represent a sampling region. It comes in three different
formats: a simple rectangle, a polygon (vector format) or a binary mask (TRUE/FALSE grid; raster format). These can all be coerced to polygonal form internally in
\pkg{spatstat} and then converted to \code{SpatialPolygons}, which is what
\code{as.SpatialPolygons.owin()} does. For objects of class
\code{SpatialPolygons} (and \code{SpatialPolygonsDataFrame}) the \pkg{sp}
polygons are simply extracted and cast into \pkg{spatstat}s polygon format
inside the \code{owin} object. For \code{SpatialPixelsDataFrame} and
\code{SpatialGridDataFrame} the grid is extracted and cast into \pkg{spatstat}s
mask format inside the \code{owin} object. In all cases any data apart from the
spatial region itself are discarded.}

\note{When coercing a SpatialPolygons object to an owin object, full 
topology checking is enabled by default. To avoid checking, set 
\code{spatstat.options(checkpolygons=FALSE)} (from spatstat (1.14-6)). 
To perform the checking later, \code{owinpolycheck(W, verbose=TRUE)}.}

\author{Edzer Pebesma \email{edzer.pebesma@uni-muenster.de}, Roger Bivand}
\examples{
run <- FALSE
if (require(spatstat.linnet, quietly=TRUE)) run <- TRUE
if (run) {
## SpatialPixelsDataFrame -> owin
data(meuse.grid) # A data.frame
gridded(meuse.grid) = ~x+y # Now a SpatialPixelsDataFrame
mg_owin <- as(meuse.grid, "owin")
mg_owin
}
if (run) {
## SpatialGridDataFrame -> owin
fullgrid(meuse.grid) <- TRUE # Now a SpatialGridDataFrame
mg_owin2 <- as(meuse.grid, "owin")
}
if (run) {
## SpatialPolygons region with a hole
ho_sp <- SpatialPolygons(list(Polygons(list(Polygon(cbind(c(0,1,1,0,0),
  c(0,0,1,1,0))), Polygon(cbind(c(0.6,0.4,0.4,0.6,0.6), 
  c(0.2,0.2,0.4,0.4,0.2)), hole=TRUE)), ID="ho")))
plot(ho_sp, col="red", pbg="pink")
}
if (run) {
## SpatialPolygons -> owin
ho <- as(ho_sp, "owin")
plot(ho)
}
if (run) {
## Define owin directly and check they are identical
ho_orig <- owin(poly=list(list(x=c(0,1,1,0), y=c(0,0,1,1)),
  list(x=c(0.6,0.4,0.4,0.6), y=c(0.2,0.2,0.4,0.4))))
identical(ho, ho_orig)
}
if (run) {
## owin -> SpatialPolygons
ho_sp1 <- as(ho, "SpatialPolygons")
all.equal(ho_sp, ho_sp1, check.attributes=FALSE)
}
}
\keyword{spatial}
