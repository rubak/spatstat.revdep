% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_result.R
\name{plot_result}
\alias{plot_result}
\title{Plot the best results}
\usage{
plot_result(
  result,
  Polygon1,
  best = 3,
  plotEn = 1,
  topographie = FALSE,
  Grid,
  Projection,
  sourceCCLRoughness,
  sourceCCL,
  weibullsrc
)
}
\arguments{
\item{result}{The output of \code{\link{windfarmGA}} or
\code{\link{genetic_algorithm}}}

\item{Polygon1}{The considered area as shapefile}

\item{best}{A numeric value indicating how many of the best individuals
should be plotted}

\item{plotEn}{A numeric value that indicates if the best energy or efficiency
output should be plotted. If (plotEn==1) plots the best energy solutions
and (plotEn==2) plots the best efficiency solutions}

\item{topographie}{A logical value, indicating whether terrain effects should
be considered and plotted or not}

\item{Grid}{The grid as SpatialPolygons, which is obtained from
\code{\link{grid_area}} and used for plotting}

\item{Projection}{A desired Projection can be used instead of the default
Lambert Azimuthal Equal Area Projection}

\item{sourceCCLRoughness}{The source to the adapted Corine Land Cover legend
as .csv file. Only required when terrain effect model is activated. As
default a .csv file within this package (\file{~/extdata/clc_legend.csv})
is taken that was already adapted manually}

\item{sourceCCL}{The source to the Corine Land Cover raster (.tif). Only
required, when the terrain effect model is activated}

\item{weibullsrc}{A list of Weibull parameter rasters, where the first list
item must be the shape parameter raster k and the second item must be the
scale parameter raster a of the Weibull distribution. If no list is given,
then rasters included in the package are used instead, which currently only
cover Austria.}
}
\value{
Returns a data.frame of the best (energy/efficiency) individual
  during all iterations
}
\description{
Plot the best solutions of the genetic algorithm.
  Depending on \code{plotEn}, either the best energy or efficiency solutions
  can be plotted. \code{best} indicates the amount of best solutions to plot.
}
\examples{
\donttest{
## Add some data examples from the package
load(file = system.file("extdata/resultrect.rda", package = "windfarmGA"))
load(file = system.file("extdata/resulthex.rda", package = "windfarmGA"))
load(file = system.file("extdata/polygon.rda", package = "windfarmGA"))

## Plot the results of a hexagonal grid optimization
Grid <- hexa_area(polygon, size = 75, FALSE)
plot_result(resulthex, polygon, best = 1, plotEn = 1, topographie = FALSE,
           Grid = Grid[[2]])

## Plot the results of a rectangular grid optimization
Grid <- grid_area(polygon, resol = 150, 1, FALSE)
plot_result(resultrect, polygon, best = 1, plotEn = 1, topographie = FALSE,
           Grid = Grid[[2]])
}
}
\seealso{
Other Plotting Functions: 
\code{\link{interpol_view}()},
\code{\link{plot_cloud}()},
\code{\link{plot_development}()},
\code{\link{plot_evolution}()},
\code{\link{plot_fitness_evolution}()},
\code{\link{plot_heatmap}()},
\code{\link{plot_parkfitness}()},
\code{\link{plot_viewshed}()},
\code{\link{plot_windfarmGA}()},
\code{\link{plot_windrose}()},
\code{\link{random_search_single}()}
}
\concept{Plotting Functions}
