\name{getKMLcoordinates}
\alias{getKMLcoordinates}

\title{Get a list of coordinates out of a KML file}

\description{
  This function parses a KML file to get the content of \code{<coordinates>} tags and returns a list of matrices representing the longitude-latitute or if \code{ignoreAltitude} is \code{FALSE} the longitude-latitute-altitude coordinates of a KML geometry.
}

\usage{
getKMLcoordinates(kmlfile, ignoreAltitude=FALSE)
}

\arguments{
  \item{kmlfile}{connection object or a character string of the KML file}
  \item{ignoreAltitude}{if set to \code{TRUE} the altitude values of a KML points will be ignored}
}

\value{
\code{coords} is a list of matrices representing the longitude-latitute or if \code{ignoreAltitude} is \code{FALSE} the longitude-latitute-altitude coordinates
}

\author{Hans-J. Bibiko}

\seealso{\code{\link{kmlPolygon}}, \code{\link{kmlLine}}}
\examples{
data(wrld_simpl)
## creates a KML file containing the polygons of South Africa (plus hole)
sw <- slot(wrld_simpl[wrld_simpl$NAME=="South Africa",], "polygons")[[1]]
tf <- tempfile()
kmlPolygon(sw, kmlfile=tf, name="South Africa", col="#df0000aa", lwd=5, 
    border=4, kmlname="R Test", 
    kmldescription="This is <b>only</b> a <a href='http://www.r-project.org'>R</a> test.")
zz <- getKMLcoordinates(tf, ignoreAltitude=TRUE)
str(zz)
zz <- getKMLcoordinates(system.file("shapes/Testing.kml", package="maptools"))
str(zz)
}
\keyword{spatial}

