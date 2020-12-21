\name{tdrankperc}
\alias{tdrankperc}
\title{percentiles for TD rank curve}
\description{gets percentiles for TD rank curve}
\usage{tdrankperc( td, percentiles = c( 0.5, 1, 2, 5, 95 ),
                   type = c( "quantile", "(i-1)/(n-1)", "i/(n+1)", "i/n" ),
                   smooth = TRUE, smoothFunction = tdrankglm )}
\arguments{
  \item{td}{table with total-deviation values}
  \item{percentiles}{percentiles at which to calculate cutoff values}
  \item{type}{see \code{\link{wtd.quantile}} for a list of different options}
  \item{smooth}{whether to use a function to smooth the results or not. Default is TRUE}
  \item{smoothFunction}{if \code{smooth} is \code{TRUE} is true, the function to use for smoothing. Default is \code{\link{tdrankglm}}, a GLM fit which was proven to do a good fit for average over subject of TD rank curves. This function is not really a smoothing procedure, but a parametric fit}
}
\details{details}

\value{
percentiles for TD rank curves
}

\seealso{ \code{\link{tdrank}}, \code{\link{tdrankglm}} }

\author{Ivan Marin-Franch}