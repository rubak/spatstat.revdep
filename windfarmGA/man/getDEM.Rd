% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visibility.R
\name{getDEM}
\alias{getDEM}
\title{Get DEM raster}
\usage{
getDEM(polygon, ISO3 = "AUT", clip = TRUE)
}
\arguments{
\item{polygon}{A Spatial / SimpleFeature Polygon to crop the DEM}

\item{ISO3}{The ISO3 code of the country}

\item{clip}{boolean, indicating if polygon should be cropped. Default is TRUE}
}
\value{
A list with the DEM raster, and a SpatialPolygonsDataFrame or NULL if
  no polygon is given
}
\description{
Get a DEM raster for a country based on ISO3 code
}
\examples{
\dontrun{
library(sp)
library(raster)
Polygon1 <- Polygon(rbind(c(4488182, 2667172), c(4488182, 2669343),
                          c(4499991, 2669343), c(4499991, 2667172)))
Polygon1 <- Polygons(list(Polygon1), 1);
Polygon1 <- SpatialPolygons(list(Polygon1))
Projection <- "+init=epsg:3035"
proj4string(Polygon1) <- CRS(Projection)
DEM_meter <- getDEM(Polygon1)
plot(DEM_meter[[1]])
plot(DEM_meter[[2]], add=TRUE)
}
}
\seealso{
Other Helper Functions: 
\code{\link{dup_coords}()},
\code{\link{getISO3}()},
\code{\link{get_grids}()},
\code{\link{grid_area}()},
\code{\link{hexa_area}()},
\code{\link{isSpatial}()},
\code{\link{permutations}()},
\code{\link{readintegerSel}()},
\code{\link{readinteger}()},
\code{\link{splitAt}()},
\code{\link{tess2SPdf}()},
\code{\link{windata_format}()}
}
\concept{Helper Functions}
