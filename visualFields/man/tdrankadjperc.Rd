\name{tdrankadjperc}
\alias{tdrankadjperc}
\title{percentiles for adjusted TD rank curve}
\description{gets percentiles for adjusted TD rank curve}
\usage{tdrankadjperc( td, percentiles = c( 0.5, 1, 2, 5, 95 ), type = "conventional",
               typequantiles = c( "quantile", "(i-1)/(n-1)", "i/(n+1)", "i/n" ),
               smooth = TRUE, smoothFunction = tdrankglm )}
\arguments{
  \item{td}{\code{\link{vf}}-object with total-deviation values}
  \item{type}{whether to use a \code{conventional} way to plot the rank TD curve or \code{ghrank} type where the \code{vf} object passed is the reconstructed within-normal TD rank curve. Default is \code{conventional}}
  \item{percentiles}{percentiles at which to calculate cutoff values}
  \item{typequantiles}{see \code{\link{wtd.quantile}} for a list of different options}
  \item{smooth}{whether to use a function to smooth the results or not. Default is TRUE}
  \item{smoothFunction}{if \code{smooth} is \code{TRUE} is true, the function to use for smoothing. Default is \code{\link{tdrankglm}}, a GLM fit which was proven to do a good fit for average over subject of TD rank curves. This function is not really a smoothing procedure, but a parametric fit}
}

\value{percentiles for adjusted TD rank curves}

\seealso{\code{\link{tdrank}}, \code{\link{tdrankperc}}}

\examples{
  \donttest{
    tdrankadjperc( td )
  }
}

\author{Ivan Marin-Franch}