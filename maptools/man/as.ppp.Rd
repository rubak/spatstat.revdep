\name{as.ppp}
\alias{as.ppp.SpatialPoints}
\alias{coerce,SpatialPoints,ppp-method}
\alias{as.ppp.SpatialPointsDataFrame}
\alias{coerce,SpatialPointsDataFrame,ppp-method}

\alias{as.SpatialPoints.ppp}
\alias{coerce,ppp,SpatialPoints-method}
\alias{as.SpatialPointsDataFrame.ppp}
\alias{coerce,ppp,SpatialPointsDataFrame-method}
\alias{as.SpatialGridDataFrame.ppp}
\alias{coerce,ppp,SpatialGridDataFrame-method}

\title{Coercion between sp objects and spatstat ppp objects}
\description{
Functions to convert between \pkg{spatstat}s planar point pattern (\code{ppp})
format and \pkg{sp}s \code{SpatialPoints} and \code{SpatialPointsDataFrame} as
well as one-way conversion from \code{SpatialGridDataFrame} to \code{ppp}.
S4-style \code{as()} coercion can be used as well.
}
\section{Methods}{
  \describe{
\item{coerce}{\code{signature(from = "SpatialPoints", to = "ppp")}}
\item{coerce}{\code{signature(from = "SpatialPointsDataFrame", to = "ppp")}}
\item{coerce}{\code{signature(from = "ppp", to = "SpatialGridDataFrame")}}
\item{coerce}{\code{signature(from = "ppp", to = "SpatialPointsDataFrame")}}
\item{coerce}{\code{signature(from = "ppp", to = "SpatialPoints")}}
}
}

\usage{
as.ppp.SpatialPoints(X)
as.ppp.SpatialPointsDataFrame(X)
as.SpatialPoints.ppp(from)
as.SpatialPointsDataFrame.ppp(from)
as.SpatialGridDataFrame.ppp(from)
}
\arguments{
  \item{from, X}{object to coerce from}
}
\details{The main conversion is between \pkg{sp}s
\code{SpatialPoints}/\code{SpatialPointsDataFrame} and \pkg{spatstat}s
\code{ppp}. Conversion between \code{SpatialGridDataFrame} and \code{ppp} should
rarely be used; using \code{\link{as.owin.SpatialGridDataFrame}} is more transparent.}

\note{The \code{ppp} format requires an observation window which is the sampling
region. The \pkg{sp} formats contain no such information and by default the
bounding box of the points is simply used. This is almost never the correct
thing to do! Rather, information about the sampling region should be converted
into \pkg{spatstat}s \code{\link[spatstat]{owin}} format and assigned as the
observation window. Usually conversion from \code{ppp} to \pkg{sp} format simply
discards the \code{owin}. However, \code{as.SpatialGridDataFrame.ppp} actually
first discards the points(!), second checks that the corresponding \code{owin}
is in a grid format (matrix of TRUE/FALSE for inside/outside sampling region),
and finally converts the TRUE/FALSE grid to a \code{SpatialGridDataFrame}.}

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

\author{Edzer Pebesma \email{edzer.pebesma@uni-muenster.de}, Roger Bivand}
\examples{
run <- FALSE
if (require(spatstat, quietly=TRUE)) run <- TRUE
if (run) {
## Convert SpatialPointsDataFrame into a marked ppp
data(meuse)
coordinates(meuse) = ~x+y
meuse_ppp <- as(meuse, "ppp")
meuse_ppp # Window is the bounding rectangle
}
if (run) {
plot(meuse_ppp, which.marks = "zinc")
}
if (run) {
## Convert SpatialPoints into an unmarked ppp
meuse2 <- as(meuse, "SpatialPoints")
as(meuse2, "ppp")
}
if (run) {
## Get sampling region in grid format and assign it as observation window
data(meuse.grid)
gridded(meuse.grid) <- ~x+y
mg_owin <- as(meuse.grid, "owin")
Window(meuse_ppp) <- mg_owin
meuse_ppp # Window is now a binary image mask (TRUE/FALSE grid)
}
if (run) {
plot(meuse_ppp, which.marks = "zinc")
}
if (run) {
## Convert marked ppp back to SpatialPointsDataFrame
rev_ppp_SPDF <- as.SpatialPointsDataFrame.ppp(meuse_ppp)
summary(rev_ppp_SPDF)
}
if (run) {
## Convert marked ppp back to SpatialPoints (discarding marks)
rev_ppp_SP <- as.SpatialPoints.ppp(meuse_ppp)
summary(rev_ppp_SP)
}
if (run) {
## Convert marked ppp back to SpatialGridDataFrame (extracting the window grid)
rev_ppp_SGDF <- as.SpatialGridDataFrame.ppp(meuse_ppp)
summary(rev_ppp_SGDF)
}
}
\keyword{spatial}
