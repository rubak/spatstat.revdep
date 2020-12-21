\name{readAsciiGrid}
\alias{readAsciiGrid}
\alias{writeAsciiGrid}
\title{ read/write to/from (ESRI) asciigrid format }
\description{ read/write to/from ESRI asciigrid format; a fuzz factor has been added to \code{writeAsciiGrid} to force cell resolution to equality if the difference is less than the square root of machine precision}
\usage{ 
readAsciiGrid(fname, as.image = FALSE, plot.image = FALSE, 
 colname = basename(fname), proj4string = CRS(as.character(NA)), 
 dec=options()$OutDec)
writeAsciiGrid(x, fname, attr = 1, na.value = -9999, dec=options()$OutDec, ...)
}

\arguments{
\item{fname}{ file name }
\item{as.image}{logical; if TRUE, a list is returned, ready to be
shown with the \code{image} command; if FALSE an object of
class \link[sp]{SpatialGridDataFrame-class} is returned }
\item{plot.image}{logical; if TRUE, an image of the map is plotted}
\item{colname}{alternative name for data column if not file basename}
\item{proj4string}{A CRS object setting the projection arguments of the Spatial Grid returned}
\item{dec}{decimal point character.  This should be a character string containing just one single-byte character --- see note below.}
\item{x}{ object of class \link[sp]{SpatialGridDataFrame}}
\item{attr}{attribute column; if missing, the first
column is taken; a name or a column number may be given}
\item{na.value}{numeric; value given to missing valued cells in the
resulting map}
\item{...}{ arguments passed to \link{write.table}, which is used
to write the numeric data }
}
\value{
	\code{readAsciiGrid} returns the grid map read; either as
	an object of class \link[sp]{SpatialGridDataFrame-class} or, if
	as.image is TRUE, as list with components \code{x}, \code{y} and \code{z}.
}
\note{In ArcGIS 8, it was not in general necessary to set the \code{dec} argument; it is not necessary in a mixed environment with ArcView 3.2 (R writes and ArcView reads "."), but inter-operation with ArcGIS 9 requires care because the defaults used by ArcGIS seem to be misleading, and it may be necessary to override what appear to be platform defaults by setting the argument.}
\seealso{ \code{\link[sp]{image}}, \code{\link[graphics]{image}} }
\author{Edzer Pebesma, edzer.pebesma@uni-muenster.de}
\examples{
x <- readAsciiGrid(system.file("grids/test.ag", package="maptools")[1])
summary(x)
image(x)
xp <- as(x, "SpatialPixelsDataFrame")
abline(h=332000, lwd=3)
xpS <- xp[coordinates(xp)[,2] < 332000,]
summary(xpS)
xS <- as(xpS, "SpatialGridDataFrame")
summary(xS)
tmpfl <- paste(tempdir(), "testS.ag", sep="/")
writeAsciiGrid(xS, tmpfl)
axS <- readAsciiGrid(tmpfl)
opar <- par(mfrow=c(1,2))
image(xS, main="before export")
image(axS, main="after import")
par(opar)
unlink(tmpfl)
}
\keyword{programming}
