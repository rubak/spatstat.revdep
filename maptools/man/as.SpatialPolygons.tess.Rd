\name{as.SpatialPolygons.tess}
\alias{as.SpatialPolygons.tess}
\alias{coerce,tess,SpatialPolygons-method}

\title{Coercion of spatstat tess object to sp SpatialPolygons object}
\description{
This function coerces \pkg{spatstat}s tessellation objects of class \code{tess}
to \pkg{sp}s \code{SpatialPolygons} class. S4-style as() coercion works as well.
}
\section{Methods}{
  \describe{
\item{coerce}{\code{signature(from = "tess", to = "SpatialPolygons")}}
}
}

\usage{
as.SpatialPolygons.tess(x)
}
\arguments{
  \item{x}{\pkg{spatstat} object of class \code{tess} to coerce from}
}

\author{Edzer Pebesma \email{edzer.pebesma@uni-muenster.de}, Roger Bivand}
\examples{
run <- FALSE
if (require(spatstat.core, quietly=TRUE)) run <- TRUE
if (run) {
A <- tess(xgrid=0:4,ygrid=0:4)
A_sp <- as(A, "SpatialPolygons")
plot(A_sp)
text(coordinates(A_sp), labels=row.names(A_sp), cex=0.6)
}
}
\keyword{spatial}
