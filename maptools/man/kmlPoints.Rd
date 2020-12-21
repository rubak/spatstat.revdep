\name{kmlPoints}
\alias{kmlPoints}

\title{Create and write a KML file on the basis of a given Points object}
\description{
  The function is used to create and write a KML file on the basis of a given SpatialPointsDataFrame object for the usage in Google Earth resp. Google Maps.
}
\usage{
kmlPoints(obj=NULL, kmlfile=NULL, kmlname="", kmldescription="",
    name=NULL, description="",
    icon="http://www.gstatic.com/mapspro/images/stock/962-wht-diamond-blank.png")
}

\arguments{
  \item{obj}{a \code{SpatialPointsDataFrame} object}
  \item{kmlfile}{if not \code{NULL} the name as character string of the kml file to be written}
  \item{kmlname}{the name of the KML layer}
  \item{kmldescription}{the description of the KML layer (HTML tags allowed)}
  \item{name}{a character vector to be used as names for each KML Placemark}
  \item{description}{a character vector to be used as the description for each KML Placemark (HTML tags allowed)}
  \item{icon}{a character vector of icon URLs to be used in the style associated with each KML Placemark}
}

\details{
The function is used to convert a given \code{SpatialPointsDataFrame} object into a series of KML Placemarks, each with a single Point. If \code{kmlfile} is not \code{NULL} the result will be written into that file. If \code{kmlfile} is \code{NULL} the generated KML lines will be returned (see also value). 

If \code{name=NULL}, the <name> tag for each Placemark will be \code{'site #'}. If a single value is used for 
\code{name} or \code{description}, that value will be replicated for each Placemark. If a single value is used for
\code{icon}, only a single style will be created and that style will be referenced by each Placemark.

Note that the geometries should be in geographical coordinates with datum WGS84.
}

\value{
x is a list with the elements \code{style} and \code{content} containing the generated lines of the KML file as character vectors if \code{kmlfile} is \code{NULL}.

y is a list with the elements \code{header} and \code{footer} representing the KML file' header resp. footer if \code{obj} is \code{NULL}.
}

\section{KML icons}{
The default icon URL is \url{http://www.gstatic.com/mapspro/images/stock/962-wht-diamond-blank.png}. Additional icons are available at:
\url{http://sites.google.com/site/gmapsdevelopment}.
}

\author{Jonathan Callahan}
\seealso{\code{\link{kmlLine}}, \code{\link{kmlOverlay}}, \code{\link{kmlPolygon}}, \code{\link[sp]{Line}}}

\examples{
data(SplashDams)
num <- length(SplashDams)
td <- tempdir()
kmlfile <- paste(td, "OregonSplashDams.kml", sep="/")
kmlname <- "Oregon Splash Dams"
kmldescription <- paste("Data for Splash Dams in western Oregon.",
 "See http://www.fs.fed.us/pnw/lwm/aem/people/burnett.html#projects_activities",
 "for more information.")
icon <- "http://www.gstatic.com/mapspro/images/stock/962-wht-diamond-blank.png"
name <- paste("Dam on",SplashDams$streamName)
description <- paste("<b>owner:</b>", SplashDams$owner, "<br><b>dates:</b>", SplashDams$datesUsed)

kmlPoints(SplashDams, kmlfile=kmlfile, name=name, description=description,
          icon=icon, kmlname=kmlname, kmldescription=kmldescription)
}
\concept{kml}
\keyword{spatial}

