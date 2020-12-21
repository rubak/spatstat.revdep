\name{LTDataFrame}
\alias{LTDataFrame}
\title{Creates covariate dataframes}
\description{
Creates covariate dataframes for observations and dummy quadrature points}
\usage{
LTDataFrame(study.area, lines, lines.psp, int.formula, det.formula, 
                   covariates, Q.lt)
}
\arguments{
  \item{study.area}{owin object that defines study area}
  \item{lines}{data frame of lines with structure as shown in \code{\link{quadscheme.lt}}}
  \item{lines.psp}{psp class with added list elements \code{width} and \code{label}}
  \item{int.formula}{model formula for intensity process}
  \item{det.formula}{model formula for detection scale process}
  \item{covariates}{covariate dataframe (see DSpat for structure)}
  \item{Q.lt}{lt quadscheme of class quad}
}
\value{ 
\item{cov.df}{dataframe of covariates followed by rows for covariates for dummy quadrature points}
\item{covariate.im}{list of covariate images}
}
\details{
Checks to make sure that all of the variables used in formula are either in \code{covariates} 
or in \code{lines}. Then it extracts the values of \code{covariates} for
each observation and for dummy points.   These are merged with the needed covariates
from \code{lines} and then a single dataframe is returned with the observations followed
by the dummy points.  In addition, the covariate images in a list are returned to keep
with the \code{dspat} object for use in \code{\link{integrate.intensity}}.
}
\author{ Devin Johnson; Jeff Laake}
\seealso{ \code{\link{quadscheme.lt}},\code{\link{create.covariate.images}} }

