\name{gloperc}
\alias{gloperc}
\title{percentiles for global indices}
\description{gets the percentiles for global indices}
\usage{gloperc( vals, percentiles = c( 0.5, 1, 2, 5, 95 ),
         type = c( "quantile", "(i-1)/(n-1)", "i/(n+1)", "i/n" ) )}
\arguments{
  \item{vals}{table with global indices for control subjects}
  \item{percentiles}{percentiles at which to calculate cutoff values}
  \item{type}{see \code{\link{wtd.quantile}} for a list of different options}
}
\details{
gets the percentiles for global indices. This function works in conjunction with \code{\link{vfstats}}. In \code{\link{vfstats}} all global indices refering to mean must start with \code{m} and all refereing to standard deviation must start with \code{s}, otherwise \code{\link{gloperc}} won't calculate the cutoff values correctly
}

\value{
percentiles for global indices
}

\seealso{\code{\link{vfindex}}, \code{\link{vfstats}}, \code{\link{vfiperc}}, \code{\link{locperc}}}

\author{Ivan Marin-Franch}