\name{getinfo.shape}
\alias{getinfo.shape}
\alias{print.shapehead}

\title{Get shapefile header information}
\description{
  Get shapefile  header information; the file should be given including its ".shp" extension, and the function will reconstruct the names of the database (dbf) file and the index (shx) file from these.
}
\usage{
getinfo.shape(filen)
\method{print}{shapehead}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filen}{name of file with *.shp extension}
  \item{x}{a shapehead list as returned by getinfo.shape}
  \item{...}{other arguments passed to print}
}
\details{
  The function calls code from shapelib to read shapefiles, a file format used by ESRI GIS software among others
}
\value{
The function returns a list of class shapehead.
}
\references{\url{http://shapelib.maptools.org/}}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}; shapelib by Frank Warmerdam}


\examples{
res <- getinfo.shape(system.file("shapes/fylk-val.shp", package="maptools")[1])
res
str(res)
}
\keyword{spatial}

