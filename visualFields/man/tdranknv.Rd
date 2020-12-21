\name{tdranknv}
\alias{tdranknv}
\title{normative values for TD rank curve}
\description{gets percentiles for TD rank curve}
\usage{tdranknv( td, smooth = TRUE, smoothFunction = tdrankglm )}
\arguments{
  \item{td}{table with total-deviation values}
  \item{smooth}{whether to use a function to smooth the results or not. Default is TRUE}
  \item{smoothFunction}{if \code{smooth} is \code{TRUE}, the function to use for smoothing. Default is \code{\link{tdrankglm}}, a GLM fit which was proven to do a good fit for average over subject of TD rank curves. This function is not really a smoothing procedure, but a parametric fit}
}
\details{details}

\value{
percentiles for TD rank curves
}

\seealso{ \code{\link{tdrank}}, \code{\link{tdrankglm}} }

\author{Ivan Marin-Franch}