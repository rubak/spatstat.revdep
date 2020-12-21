\name{ContourLines2SLDF}
\alias{ArcObj2SLDF}
\alias{ContourLines2SLDF}
\alias{MapGen2SL}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Converter functions to build SpatialLinesDataFrame objects}
\description{
  These functions show how to build converters to SpatialLinesDataFrame objects: 
\code{ArcObj2SLDF} from the list returned by the \code{get.arcdata} function in the RArcInfo package; \code{ContourLines2SLDF} from the list returned by the \code{contourLines} function in the graphics package (here the data frame is just the contour levels, with one Lines object made up of at least one Line object per level); and \code{MapGen2SL} reads a file in "Mapgen" format into a \code{SpatialLines} object.
}
\usage{
ArcObj2SLDF(arc, proj4string=CRS(as.character(NA)), IDs)
ContourLines2SLDF(cL, proj4string=CRS(as.character(NA)))
MapGen2SL(file, proj4string=CRS(as.character(NA)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{arc}{a list returned by the \code{get.arcdata} function in the RArcInfo package}
  \item{IDs}{vector of unique character identifiers; if not given, suitable defaults will be used, and the same values inserted as data slot row names}
  \item{cL}{a list returned by the \code{contourLines} function in the graphics package}
  \item{proj4string}{Object of class \code{"CRS"}; see \link[sp]{CRS-class}}
  \item{file}{filename of a file containing a Mapgen line data set}
}

\value{
  A SpatialLinesDataFrame object
}

\author{ Roger Bivand; Edzer Pebesma }

\note{Coastlines of varying resolution may be chosen online and downloaded in "Mapgen" text format from \url{https://www.ngdc.noaa.gov/mgg/shorelines/shorelines.html}, most conveniently using the interactive selection tool, but please note the 500,000 point limit on downloads, which is easy to exceed.}


\seealso{ \link[sp]{SpatialLines-class} }
\examples{
#data(co37_d90_arc) # retrieved as: 
# library(RArcInfo)
# fl <- "http://www.census.gov/geo/cob/bdy/co/co90e00/co37_d90_e00.zip"
# download.file(fl, "co37_d90_e00.zip")
# e00 <- zip.file.extract("co37_d90.e00", "co37_d90_e00.zip")
# e00toavc(e00, "ncar")
# arc <- get.arcdata(".", "ncar")
#res <- arcobj2SLDF(arc)
#plot(res)
#invisible(title(""))
res <- ContourLines2SLDF(contourLines(volcano))
plot(res, col=terrain.colors(nrow(as(res, "data.frame"))))
title("Volcano contours as SpatialLines")
}
\keyword{spatial}
