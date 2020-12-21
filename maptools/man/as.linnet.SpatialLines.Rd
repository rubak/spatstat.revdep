\name{as.linnet.SpatialLines}
\alias{as.linnet.SpatialLines}
\alias{coerce,SpatialLines,linnet-method}
\alias{coerce,SpatialLinesDataFrame,linnet-method}

\title{
  Convert SpatialLines to Linear Network
}
\description{
  Convert an object of class \code{SpatialLines} or
  \code{SpatialLinesDataFrame} (from package \pkg{sp}),
  representing a collection of polygonal lines,
  into an object of class \code{linnet} (from package \pkg{spatstat}),
  representing a network.
}
\usage{
  as.linnet.SpatialLines(X, \dots, fuse = TRUE)
  \S4method{coerce}{SpatialLines,linnet}(from, to = "linnet", strict = TRUE)
  \S4method{coerce}{SpatialLinesDataFrame,linnet}(from, to = "linnet",
    strict = TRUE)
}
\arguments{
  \item{X, from}{
    Object of class \code{SpatialLines} or
    \code{SpatialLinesDataFrame} to be converted.
  }
  \item{to}{output object of class \dQuote{linnet}.}
  \item{strict}{logical flag. If TRUE, the returned object must be strictly from the target class.}
  \item{\dots}{
    Ignored.
  }
  \item{fuse}{
    Logical value indicating whether to join different curves
    which have a common vertex.
  }
}
\details{
  This function converts an object of class \code{SpatialLines}
  or \code{SpatialLinesDataFrame} into an object of class \code{linnet}.
  It is not a method for the \pkg{spatstat} generic function
  \code{\link[spatstat]{as.linnet}}, but like other S4 coercion functions for
  \pkg{sp} classes to \pkg{spatstat} classes, it may be called directly as a function.

  An object of class {SpatialLines}
  or \code{SpatialLinesDataFrame} (from package \pkg{sp}) represents
  a list of lists of the coordinates of lines, such as a list of all roads in
  a city. An object of class \code{linnet} in the \pkg{spatstat} package
  represents a linear network, such as a road network.

  If \code{fuse=FALSE}, each \dQuote{Line} object in \code{X}
  will be treated as if it were disconnected from the others.
  The result is a network that consists of many disconnected sub-networks,
  equivalent to the list of \dQuote{Line} objects.
  
  If \code{fuse=TRUE} (the default), the code will
  search for identical pairs of vertices occurring in different
  \dQuote{Line} objects, and will treat them as identical vertices,
  effectively joining the two \dQuote{Line} objects at this common
  vertex.

  If \code{X} belongs to class \code{SpatialLinesDataFrame},
  the associated columns of data in the auxiliary data frame
  \code{slot(X, "data")} will be copied to the output as the marks
  attached to the line segments of the network. See the Examples.
}
\value{
  An object of class \code{linnet}.
}
\author{
  Adrian Baddeley.
}
\seealso{
  \code{\link[spatstat]{as.linnet}}
}
\examples{
  run <- FALSE
  if(require(spatstat, quietly=TRUE)) run <- TRUE
if (run) {
   dname <- system.file("shapes", package="maptools")
   fname <- file.path(dname, "fylk-val.shp")
   fylk <- readShapeSpatial(fname, proj4string=CRS("+proj=utm +zone=33 +ellps=WGS84"))
   is.projected(fylk)
}
if (run) {
   L <- as(fylk, "linnet")
   print(max(vertexdegree(L)))
}
if (run) {
   L0 <- as.linnet.SpatialLines(fylk, fuse=FALSE)
   print(max(vertexdegree(L0)))
}
if (run) {
   ## extract data associated with each network segment
   head(marks(as.psp(L)))
}
if (run) {
   fname <- file.path(dname, "fylk-val-ll.shp")
   fylk <- readShapeSpatial(fname, proj4string=CRS("+proj=longlat +ellps=WGS84"))
   is.projected(fylk)
}
if (run) {
   try(L <- as(fylk, "linnet"))
  }
}
\keyword{spatial}
\keyword{manip}
