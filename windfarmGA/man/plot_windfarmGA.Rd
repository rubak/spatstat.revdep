% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_windfarmGA.R
\name{plot_windfarmGA}
\alias{plot_windfarmGA}
\title{Plot the results of an optimization run}
\usage{
plot_windfarmGA(
  result,
  Polygon1,
  GridMethod = "r",
  whichPl = "all",
  best = 1,
  plotEn = 1,
  Projection,
  weibullsrc
)
}
\arguments{
\item{result}{The output of \code{\link{windfarmGA}} or
\code{\link{genetic_algorithm}}}

\item{Polygon1}{The area as shapefile.}

\item{GridMethod}{Which grid spacing method was used. Default is
"rectangular". If hexagonal grid cells were used, assign any of the
following arguments: "h","hexa", "hexagonal"}

\item{whichPl}{Which plots should be shown: 1-6 are possible. The default is
"all" which shows all available plots}

\item{best}{A numeric value indicating how many of the best individuals
should be plotted}

\item{plotEn}{A numeric value that indicates if the best energy or efficiency
output is plotted. If (plotEn==1) plots the best energy solutions
and (plotEn==2) plots the best efficiency solutions}

\item{Projection}{A desired Projection can be used instead of the default
Lambert Azimuthal Equal Area Projection}

\item{weibullsrc}{A list of Weibull parameter rasters, where the first list
item must be the shape parameter raster `k` and the second item must be the
scale parameter raster `a` of the Weibull distribution. If no list is given,
then rasters included in the package are used instead, which currently only
cover Austria.}
}
\description{
Plot the results of a genetic algorithm run with given inputs.
  Several plots try to show all relevant effects and outcomes of the
  algorithm. 6 plot methods are available that can be selected individually.
}
\examples{
\donttest{
library(sp)
## Add some data examples from the package
load(file = system.file("extdata/resultrect.rda", package = "windfarmGA"))
load(file = system.file("extdata/resulthex.rda", package = "windfarmGA"))
load(file = system.file("extdata/polygon.rda", package = "windfarmGA"))

## Plot the results of a hexagonal grid optimization
plot_windfarmGA(resulthex, GridMethod = "h", polygon, whichPl = "all", best = 1, plotEn = 1)

## Plot the results of a rectangular grid optimization
plot_windfarmGA(resultrect, GridMethod = "r", polygon, whichPl = "all", best = 1, plotEn = 1)
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
\code{\link{plot_result}()},
\code{\link{plot_viewshed}()},
\code{\link{plot_windrose}()},
\code{\link{random_search_single}()}
}
\concept{Plotting Functions}
