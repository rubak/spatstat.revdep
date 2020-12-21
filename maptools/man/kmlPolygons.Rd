\name{kmlPolygons}
\alias{kmlPolygons}

\title{Create and write a KML file on the basis of a given Polygons object or list of Polygons or SpatialPolygonsDataFrame }

\description{
The function is used to create and write a KML file on the basis of a given Polygons object (a list of Polygon objects of SpatialPolygonsDataFrame class) for the usage in Google Earth and Google Maps.
}

\usage{
kmlPolygons(obj=NULL, kmlfile=NULL, 
    name="KML Polygons", description="", col=NULL, visibility=1, lwd=1,
    border="white", kmlname="", kmldescription="")
}

\arguments{
  \item{obj}{a \code{Polygons} or \code{SpatialPolygonsDataFrame} object or list of objects}
  \item{kmlfile}{if not \code{NULL} the name as character string of the kml file to be written to working directory as "NAME.kml"}
  \item{name}{the name of the KML polygon in Google Earth}
  \item{description}{the description of the KML polygon displayed in Google Earth or Maps (HTML tags allowed)}
  \item{col}{the fill color (see also Color Specification) of the KML polygon. If passing a list of \code{Polyons} or \code{SpatialPolygonsDataFrame} and \code{length(col)} is less than \code{length(object)} the first color in col will be applied to all objects in the list}
  \item{visibility}{if set to \code{1} or \code{TRUE} specifies that the KML polygon should be visible after loading}
  \item{lwd}{the stroke (polygon's border line) width for the KML polygon}
  \item{border}{the stroke color (see also Color Specification) for the KML polygon}
  \item{kmlname}{the name of the KML layer}
  \item{kmldescription}{the description of the KML layer (HTML tags allowed)}
}

\details{
    The function is used to convert a given \code{Polygons} object (a list of Polygon objects) or the  \code{Polygons} object listed in a passed \code{SpatialPolygonsDataFrame} object into KML polygon. If \code{kmlfile} is not \code{NULL} the result will be written into that file. If \code{kmlfile} is \code{NULL} the generated KML lines will be returned (see also value).

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

\author{Hans-J. Bibiko, Jon Callihan, Steven Brey}

\seealso{\code{\link{kmlPolygon}}, \code{\link{kmlLines}}, \code{\link[sp]{SpatialPolygons}}, \code{kmlPoints}}

\examples{
data(wrld_simpl)
td <- tempdir()
kmlfile <- paste(td, "worldPolitical.kml", sep="/")
## creates a KML file containing the polygons of a political world map
kmlPolygons(wrld_simpl, kmlfile = kmlfile, name = "KML Polygons", 
         description = "the world", col = "red",
         visibility = 1, lwd = 1, border = "white", kmlname = "R Test", 
         kmldescription = "This is <b>only</b> a <a href='http://www.r-project.org'>R</a> test.")

data(wrld_simpl)
## create a KML file containing the polygons of Brazil, Uganda, and Canada
regions <- c("Brazil","Canada","Uganda")
wrld_simpl_subset <- wrld_simpl[wrld_simpl$NAME \%in\% regions,]
kmlfile <- paste(td, "worldPoliticalSubset.kml", sep="/")
kmlPolygons(wrld_simpl_subset, kmlfile = kmlfile,
 name = "KML Polygons subset", description = "three countries", col = "blue", 
 visibility = 1, lwd = 1, border = "white", kmlname = "R Test 2", 
 kmldescription = "This is <b>only</b> a <a href='http://www.r-project.org'>R</a> test.")
## combine to make a list of polygon objects to plot
polList <- c(regions,wrld_simpl)
kmlfile <- paste(td, "worldPoliticalandSubset.kml", sep="/")
kmlPolygons(wrld_simpl_subset, kmlfile = kmlfile,
 name = "KML Polygons subset", description = "three countries highlighted in world", 
 col = sample(colours(), length(polList)), visibility = 1, lwd = 1, border = "white", 
 kmlname = "R Test 2", 
 kmldescription = "This is <b>only</b> a <a href='http://www.r-project.org'>R</a> test.")

}

\concept{kml}
\keyword{spatial}

