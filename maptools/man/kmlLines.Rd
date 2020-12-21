\name{kmlLines}
\alias{kmlLines}

\title{Create and write a KML file on the basis of a given Lines object}
\description{
  The function is used to create and write a KML file on the basis of a given Lines object (a list of Line objects) for the usage in Google Earth and Google Maps.
}
\usage{
kmlLines(obj=NULL, kmlfile=NULL, 
        name="R Lines", description="", col=NULL, visibility=1, lwd=1,
        kmlname="", kmldescription="")
}

\arguments{
  \item{obj}{a \code{Lines} or \code{SpatialLinesDataFrame} object}
  \item{kmlfile}{if not \code{NULL} the name as character string of the kml file to be written}
  \item{name}{the name of the KML line}
  \item{description}{the description of the KML line (HTML tags allowed)}
  \item{col}{the stroke color (see also Color Specification) of the KML line}
  \item{visibility}{if set to \code{1} or \code{TRUE} specifies that the KML line should be visible after loading}
  \item{lwd}{the stroke width for the KML line}
  \item{kmlname}{the name of the KML layer}
  \item{kmldescription}{the description of the KML layer (HTML tags allowed)}
}

\details{
  The function is used to convert a given \code{Lines} object (a list of Line objects) or the first \code{Lines} object listed in a passed \code{SpatialLinesDataFrame} object into KML line(s). If \code{kmlfile} is not \code{NULL} the result will be written into that file. If \code{kmlfile} is \code{NULL} the generated KML lines will be returned (see also value). Function no longer uses \code{append} greatly improving performance on large objects or lists. 

For a passed \code{Lines} object the function generates a <Style> tag whereby its id attribute is set to the passed object's ID.

Note that the geometries should be in geographical coordinates with datum WGS84.

The resulting KML line will be embedded in \code{<Placemark><MultiGeometry><LineString>}.
}

\value{
x is a list with the elements \code{style} and \code{content} containing the generated lines of the KML file as character vectors if \code{kmlfile} is \code{NULL}.

y is a list with the elements \code{header} and \code{footer} representing the KML file header and footer if \code{obj} is \code{NULL}.
}

\section{Color Specification}{
The following color specifications are allowed: \code{'red'}, \code{2}, or as hex code \code{'#RRGGBB'} resp. \code{'#RRGGBBAA'} for passing the alpha value.
}

\author{Hans-J. Bibiko, Jon Callahan, Steven Brey}
\seealso{\code{\link{kmlOverlay}}, \code{\link{kmlPolygon}}, \code{\link[sp]{Line}}}

\examples{
# Maptools library required
library(maptools)
# load line object
rivers <- readShapeSpatial(system.file("shapes/fylk-val-ll.shp",
                           package="maptools")[1], proj4string=CRS("+proj=longlat +ellps=WGS84"))
# create kml file
td <- tempdir()
kmlfile <- paste(td, "rivers.kml", sep="/")
kmlLines(rivers, kmlfile = kmlfile, name = "R Lines", 
         description = "Hello!", col = "blue", visibility = 1, lwd = 1, 
         kmlname = "", kmldescription = "")
}
\concept{kml}
\keyword{spatial}

