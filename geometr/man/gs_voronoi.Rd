% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gs_voronoi.R
\name{gs_voronoi}
\alias{gs_voronoi}
\title{Create a voronoi tiling \code{geom}}
\usage{
gs_voronoi(anchor = NULL, window = NULL, features = 3, sketch = NULL, ...)
}
\arguments{
\item{anchor}{[\code{geom(1)}|\code{data.frame(1)}]\cr Object to derive the
\code{geom} from. It must include column names \code{x}, \code{y} and
optionally a custom \code{fid}.}

\item{window}{[\code{data.frame(1)}]\cr in case the reference window deviates
from the bounding box of \code{anchor} (minimum and maximum values),
specify this here.}

\item{features}{[\code{integerish(1)}]\cr number of tiles to sketch.}

\item{sketch}{[\code{RasterLayer(1)} | \code{matrix(1)}]\cr Gridded object
that serves as template to sketch the tiling.}

\item{...}{[various]\cr graphical parameters to \code{\link{gt_locate}}, in
case the tiling is sketched; see \code{\link{gpar}}.}
}
\value{
An invisible \code{geom}.
}
\description{
Create a voronoi tiling \code{geom}
}
\examples{
# 1. create voronoi polygons programmatically
coords <- data.frame(x = c(40, 70, 70, 50),
                     y = c(40, 40, 60, 70))
window <- data.frame(x = c(0, 80),
                     y = c(0, 80))
aGeom <- gs_point(anchor = coords, window = window)
visualise(voronoi = aGeom)

tiles <- gs_voronoi(anchor = aGeom)
visualise(tiles, new = FALSE)
\donttest{
# 2. sketch a voronoi polygon by clicking into a template
gs_voronoi(sketch = gtRasters$continuous) \%>\%
  visualise(tiles = ., new = FALSE)
}
}
\seealso{
Other tilings: 
\code{\link{gs_tiles}()}
}
\concept{tilings}
