% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/windfarmGA.R
\name{windfarmGA}
\alias{windfarmGA}
\title{Run a Genetic Algorithm to optimize a wind farm layout}
\usage{
windfarmGA(
  dns,
  layer,
  Polygon1,
  GridMethod,
  Projection,
  sourceCCL,
  sourceCCLRoughness,
  vdirspe,
  Rotor = 30,
  fcrR = 3,
  n = 10,
  topograp = FALSE,
  iteration = 20,
  referenceHeight = 50,
  RotorHeight = 50,
  SurfaceRoughness = 0.14,
  Proportionality = 1,
  mutr = 0.008,
  elitism = TRUE,
  nelit = 7,
  selstate = "FIX",
  crossPart1 = "EQU",
  trimForce = TRUE,
  weibull,
  weibullsrc,
  Parallel,
  numCluster,
  verbose = FALSE,
  plotit = FALSE
)
}
\arguments{
\item{dns}{The data source name (interpretation varies by driver — for some
drivers, dsn is a file name, but may also be a folder)}

\item{layer}{The layer name}

\item{Polygon1}{The considered area as SpatialPolygon, SimpleFeature Polygon
or coordinates as matrix/data.frame}

\item{GridMethod}{Should the polygon be divided into rectangular or hexagonal
grid cells? The default is "Rectangular" grid cells and hexagonal grid
cells are computed when assigning "h" or "hexagon" to this input variable.}

\item{Projection}{A desired Projection can be used instead of the default
Lambert Azimuthal Equal Area Projection (EPSG:3035).}

\item{sourceCCL}{The path to the Corine Land Cover raster (.tif). Only
required when the terrain effect model is activated. If nothing is assign,
it will try to download a version from the EEA-website.}

\item{sourceCCLRoughness}{The source to the adapted Corine Land Cover legend
as .csv file. Only required when terrain effect model is activated. As
default a .csv file within this package (\file{~/extdata}) is taken that
was already adapted manually. To use your own .csv legend this variable has
to be assigned.}

\item{vdirspe}{A data.frame containing the incoming wind speeds, wind
directions and probabilities}

\item{Rotor}{A numeric value that gives the rotor radius in meter}

\item{fcrR}{A numeric value that is used for grid spacing. Default is 5}

\item{n}{A numeric value indicating the required amount of turbines}

\item{topograp}{Logical value, which indicates if the terrain effect model
should be enabled or not. Default is FALSE}

\item{iteration}{A numeric value indicating the desired amount of iterations
of the algorithm. Default is 20}

\item{referenceHeight}{The height at which the incoming wind speeds were
measured. Default is the RotorHeight.}

\item{RotorHeight}{The desired height of the turbine.}

\item{SurfaceRoughness}{A surface roughness length of the considered area in
m.  If the terrain effect model is activated, a surface roughness will be
calculated for every grid cell with the elevation and land cover
information. Default is 0.3}

\item{Proportionality}{A numeric value used for grid calculation. Determines
the percentage a grid has to overlay. Default is 1}

\item{mutr}{A numeric mutation rate with a default value of 0.008}

\item{elitism}{Boolean value, which indicates whether elitism should be
activated or not. Default is TRUE}

\item{nelit}{If \code{elitism} is TRUE, this input determines the amount 
of individuals in the elite group. Default is 7}

\item{selstate}{Determines which selection method is used, "FIX" selects a
constant percentage and "VAR" selects a variable percentage, depending on
the development of the fitness values. Default is "FIX"}

\item{crossPart1}{Determines which crossover method is used, "EQU" divides
the genetic code at equal intervals and "RAN" divides the genetic code at
random locations. Default is "EQU"}

\item{trimForce}{If activated (\code{trimForce == TRUE}), the algorithm will
take a probabilistic approach to trim the windfarms to the desired amount
of turbines. If deactivated (\code{trimForce == FALSE}) the adjustment will
be random. Default is FALSE}

\item{weibull}{A logical value that specifies whether to take Weibull
parameters into account. If `weibull == TRUE`, the wind speed values from
the `vdirspe` data frame are ignored. The algorithm will calculate the mean
wind speed for every wind turbine according to the Weibull parameters.
Default is FALSE}

\item{weibullsrc}{A list of Weibull parameter rasters, where the first list
item must be the shape parameter raster `k` and the second item must be the
scale parameter raster `a` of the Weibull distribution. If no list is given,
then rasters included in the package are used instead, which currently only
cover Austria. This variable is only used if `weibull == TRUE`.}

\item{Parallel}{Boolean value, indicating whether parallel processing should
be used. The parallel and doParallel packages are used for parallel
processing. Default is FALSE}

\item{numCluster}{If \code{Parallel} is TRUE, this variable defines the 
number of clusters to be used}

\item{verbose}{If TRUE it will print information for every generation.
Default is FALSE}

\item{plotit}{If TRUE it will plot the best windfarm of every generation. 
Default is FALSE}
}
\value{
The result is a matrix with aggregated values per generation, the
  best individual regarding energy and efficiency per generation, some fuzzy
  control variables per generation, a list of all fitness values per
  generation, the amount of individuals after each process, a matrix of all
  energy, efficiency and fitness values per generation, the selection and
  crossover parameters, a matrix with the generational difference in maximum
  and mean energy output, a matrix with the given inputs, a dataframe with
  the wind information, the mutation rate per generation and a matrix with all
  tested wind farm layouts.
}
\description{
The initiating function of an optimization run which will
  interactively check user-inputs. If all inputs are correct, an optimization
  will be started.
}
\details{
A terrain effect model can be included in the optimization process.
  Therefore, an SRTM elevation model will be downloaded automatically via the
  \code{raster::getData} function. A land cover raster can also be downloaded
  automatically from the EEA-website, or the path to a raster file can be
  passed to \code{sourceCCL}. The algorithm uses an adapted version of the
  Raster legend ("clc_legend.csv"), which is stored in the package directory
  \file{~/inst/extdata}. To use other values for the land cover roughness
  lengths, insert a column named \strong{"Rauhigkeit_z"} to the .csv file,
  assign a surface roughness length to all land cover types. Be sure that all
  rows are filled with numeric values and save the file with \strong{";"}
  separation. Assign the path of the file to the input variable
  \code{sourceCCLRoughness} of this function.
}
\seealso{
\code{\link{genetic_algorithm}}

Other Genetic Algorithm Functions: 
\code{\link{crossover}()},
\code{\link{fitness}()},
\code{\link{genetic_algorithm}()},
\code{\link{init_population}()},
\code{\link{mutation}()},
\code{\link{selection}()},
\code{\link{trimton}()}
}
\concept{Genetic Algorithm Functions}
