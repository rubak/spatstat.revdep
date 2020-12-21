\name{sp2tmap}
\alias{sp2tmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert SpatialPolygons object for Stata tmap command}
\description{
  The function converts a SpatialPolygons object for use with the Stata tmap command, by creating a data frame with the required columns.
}
\usage{
sp2tmap(SP)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SP}{a SpatialPolygons object}
}

\value{
 a data frame with three columns:
  \item{\_ID}{an integer vector of polygon identifiers in numeric order}
  \item{\_X}{numeric x coordinate}
  \item{\_Y}{numeric y coordinate}
  and an \code{ID_n} attribute with the named polygon identifiers
}
\references{\url{https://www.stata.com/search.cgi?query=tmap}}
\author{Roger Bivand}
\seealso{\code{\link[foreign]{write.dta}}}
\examples{
\dontrun{
xx <- readShapePoly(system.file("shapes/sids.shp", package="maptools")[1],
 IDvar="FIPSNO", proj4string=CRS("+proj=longlat +ellps=clrk66"))
plot(xx, border="blue", axes=TRUE, las=1)
tmapdf <- sp2tmap(as(xx, "SpatialPolygons"))
if (require(foreign)) {
td <- tempdir()
write.dta(tmapdf, file=file.path(td, "NCmap.dta"), version=7)
NCdf <- as(xx, "data.frame")
NCdf$ID_n <- attr(tmapdf, "ID_names")
write.dta(NCdf, file=file.path(td, "NC.dta"), version=7)
}
}
}
\keyword{spatial}

