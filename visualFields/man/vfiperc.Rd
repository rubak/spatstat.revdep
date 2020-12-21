\name{vfiperc}
\alias{vfiperc}
\title{percentiles for visual field index}
\description{gets the percentiles for visual field index}
\usage{vfiperc( vals, percentiles = c( 0.5, 1, 2, 5, 95 ),
         type = c( "quantile", "(i-1)/(n-1)", "i/(n+1)", "i/n" ) )}
\arguments{
  \item{vals}{\code{\link{vf}}-object with vfi values}
  \item{percentiles}{percentiles at which to calculate cutoff values}
  \item{type}{see \code{\link{wtd.quantile}} for a list of different options}
}

\value{percentiles for visual field index}

\references{
[1] B. Bengtsson and A. Heijl. \emph{A visual field index for calculation of glaucoma rate of progression}. American Journal of Ophthalmology, 145, 2008
}

\seealso{ \code{\link{vfindex}}, \code{\link{vfstats}}, \code{\link{gloperc}}, \code{\link{locperc}} }

\author{Ivan Marin-Franch}