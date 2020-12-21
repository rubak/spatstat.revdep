\name{thinnedSpatialPoly}
\alias{thinnedSpatialPoly}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Douglas-Peuker line generalization for Spatial Polygons}

\description{
The function applies the implementation of the Douglas-Peuker algorithm for line generalization or simplification (originally from shapefiles) to objects inheriting from Spatial Polygons. It does not preserve topology, so is suitable for visualisation, but not for the subsequent analysis of the polygon boundaries, as artefacts may be created, and boundaries of neighbouring entities may be generalized differently. If the rgeos package is available, \code{thinnedSpatialPolyGEOS} will be used with partial topology preservation instead of the R implementation here by passing arguments through.
}
\usage{
thinnedSpatialPoly(SP, tolerance, minarea=0, topologyPreserve = FALSE,
                 avoidGEOS = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SP}{an object inheriting from the SpatialPolygons class}
  \item{tolerance}{the tolerance value in the metric of the input object}
  \item{minarea}{the smallest area of Polygon objects to be retained, ignored if \pkg{rgeos} used}
  \item{topologyPreserve}{choose between two \pkg{rgeos} options: logical determining if the algorithm should attempt to preserve the topology (nodes not complete edges) of the original geometry}
  \item{avoidGEOS}{use R DP code even if \pkg{rgeos} is available}
}
\value{An object of the same class as the input object}
\references{Douglas, D. and Peucker, T. (1973). Algorithms for 
the reduction of the number of points required to 
represent a digitized line or its caricature. 
\emph{The Canadian Cartographer} 10(2). 112-122.}
\author{Ben Stabler, Michael Friendly, Roger Bivand}

\note{Warnings reporting: Non-finite label point detected and replaced, reflect the changes in the geometries of the polygons induced by line generalization.}

%\seealso{\code{\link[rgeos]{gSimplify}}}

\examples{
xx <- readShapeSpatial(system.file("shapes/sids.shp", package="maptools")[1],
      IDvar="FIPSNO", proj4string=CRS("+proj=longlat +ellps=clrk66"))
object.size(as(xx, "SpatialPolygons"))
xxx <- thinnedSpatialPoly(xx, tolerance=0.05, minarea=0.001)
object.size(as(xxx, "SpatialPolygons"))
par(mfrow=c(2,1))
plot(xx)
plot(xxx)
par(mfrow=c(1,1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

