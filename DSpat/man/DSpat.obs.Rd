\name{DSpat.obs}
\alias{DSpat.obs}
\docType{data}
\title{Observation dataframe for DSpat}
\description{
An example dataframe for the observations for fitting a spatial model with DSpat}
\usage{data(DSpat.obs)}
\format{
  A data frame with 395 observations on the following 6 variables.
  \describe{
    \item{\code{label}}{unique line label that links points to lines}
    \item{\code{x}}{x coordinate of observation point}
    \item{\code{y}}{y coordinate of observation point}
  }
}
\details{
Example observation dataframe simulated with specific covariates across a 100x100 window.
Only the x,y coordinates and line label are needed.  The covariates are
extracted based on the x,y coordinates.
}
\keyword{datasets}
