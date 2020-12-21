\name{kmlPolygon}
\alias{kmlPolygon}

\title{Create and write a KML file on the basis of a given Polygons object}

\description{
 The function is used to create and write a KML file on the basis of a given Polygons object (a list of Polygon objects) for the usage in Google Earth resp. Google Maps.
}

\usage{
kmlPolygon(obj=NULL, kmlfile=NULL, 
    name="R Polygon", description="", col=NULL, visibility=1, lwd=1, border=1,
    kmlname="", kmldescription="")
}

\arguments{
  \item{obj}{a \code{Polygons} or \code{SpatialPolygonsDataFrame} object}
  \item{kmlfile}{if not \code{NULL} the name as character string of the kml file to be written}
  \item{name}{the name of the KML polygon}
  \item{description}{the description of the KML polygon (HTML tags allowed)}
  \item{col}{the fill color (see also Color Specification) of the KML polygon}
  \item{visibility}{if set to \code{1} or \code{TRUE} specifies that the KML polygon should be visible after loading}
  \item{lwd}{the stroke width for the KML polygon}
  \item{border}{the stroke color (see also Color Specification) for the KML polygon}
  \item{kmlname}{the name of the KML layer}
  \item{kmldescription}{the description of the KML layer (HTML tags allowed)}
}

\details{
    The function is used to convert a given \code{Polygons} object (a list of Polygon objects) or the first \code{Polygons} object listed in a passed \code{SpatialPolygonsDataFrame} object into KML polygon. If \code{kmlfile} is not \code{NULL} the result will be written into that file. If \code{kmlfile} is \code{NULL} the generated KML lines will be returned (see also value). 

The conversion can also handle polygons which are marked as holes inside of the Polygons object if these holes are listed right after that polygon in which these holes appear. That implies that a given plot order set in the Polygons object will \bold{not} be considered.

For a passed \code{Polygons} object the function generates a <Style> tag whereby its id attribute is set to the passed object's ID.

Note that the geometries should be in geographical coordinates with datum WGS84.

The resulting KML polygon will be embedded in \code{<Placemark><MultiGeometry><Polygon>}.
}

\value{
x is a list with the elements \code{style} and \code{content} containing the generated lines of the KML file as character vectors if \code{kmlfile} is \code{NULL}.

y is a list with the elements \code{header} and \code{footer} representing the KML file' header resp. footer if \code{obj} is \code{NULL} (see second example).
}

\section{Color Specification}{
The following color specifications are allowed: \code{'red'}, \code{2}, or as hex code \code{'#RRGGBB'} resp. \code{'#RRGGBBAA'} for passing the alpha value.
}

\author{Hans-J. Bibiko}

\seealso{\code{\link{kmlOverlay}}, \code{\link{kmlLine}}, \code{\link[sp]{SpatialPolygons}}}

\examples{
data(wrld_simpl)
## creates a KML file containing the polygons of South Africa (plus hole)
sw <- slot(wrld_simpl[wrld_simpl$NAME=="South Africa",], "polygons")[[1]]
tf <- tempfile()
kmlPolygon(sw, kmlfile=tf, name="South Africa", col="#df0000aa", lwd=5, 
    border=4, kmlname="R Test", 
    kmldescription="This is <b>only</b> a <a href='http://www.r-project.org'>R</a> test.")
tf

## creates a KML file containing the polygons of South Africa, Switzerland, and Canada
sw  <- wrld_simpl[wrld_simpl$NAME \%in\% c("South Africa", "Switzerland", "Canada"),]
out <- sapply(slot(sw, "polygons"), function(x) { kmlPolygon(x,
    name=as(sw, "data.frame")[slot(x, "ID"), "NAME"], 
    col="red", lwd=1.5, border='black', 
    description=paste("ISO3:", slot(x, "ID"))) })
tf <- tempfile()
kmlFile <- file(tf, "w")
tf
cat(kmlPolygon(kmlname="R Test", kmldescription="<i>Hello</i>")$header, 
    file=kmlFile, sep="\n")
cat(unlist(out["style",]), file=kmlFile, sep="\n")
cat(unlist(out["content",]), file=kmlFile, sep="\n")
cat(kmlPolygon()$footer, file=kmlFile, sep="\n")
close(kmlFile)
}

\concept{kml}
\keyword{spatial}

