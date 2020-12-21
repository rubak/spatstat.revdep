\name{DSpat.covariates}
\alias{DSpat.covariates}
\docType{data}
\title{Raster covariates study area}
\description{
Example set of raster covariates for computing predicted intensity/abundance across
the entire study area (100x100).
}
\usage{data(DSpat.covariates)}
\format{
  A data frame with 10000 (1x1 raster element) on the following 4 variables.
  \describe{
    \item{\code{x}}{x coordinate for mid-point of raster cell}
    \item{\code{y}}{y coordinate for mid-point of raster cell}
    \item{\code{river}}{distance from river to center of raster element}
    \item{\code{habitat}}{type of habitat for raster element}
  }
}
\keyword{datasets}
