\name{locperc}
\alias{locperc}
\title{percentiles for each location}
\description{gets the percentiles for each location}
\usage{locperc( vals, stds, percentiles = c( 0.5, 1, 2, 5, 95 ),
         type = c( "quantile", "(i-1)/(n-1)", "i/(n+1)", "i/n" ),
         poolLocations = FALSE )}
\arguments{
  \item{vals}{table with sensitivity values, total-deviation values, or pattern-deviation values}
  \item{stds}{standard deviations per location for sensitivity values, total-deviation values, or pattern-deviation values}
  \item{percentiles}{percentiles at which to calculate cutoff values}
  \item{type}{see \code{\link{wtd.quantile}} for a list of different options}
  \item{poolLocations}{assume that the shape of the empirical distributions at each location is the same and pool sensitivity, TD, or PD values per location. This is useful when the number of controls available is small. Sample size is increased at the expense of posible bias due to the fact that empirical distributions are not necessarily the same in each location}
}
\details{
gets the percentiles at each location. Since the number of visits per subject can be variable. To account for that, weighted quantile is used in which the weights for is the inverse of the number of visits for the subject

All the data passed to the function must belong to the same perimeter \code{tperimetry}, the same pattern of locations \code{talgorithm}, and the same presentation algorithm \code{tpattern}.
}

\value{
percentiles for each location
}

\references{
[1] A. Heijl, G. Lindgren, and J. Olsson. \emph{A package for the statistical analysis of visual fields}. Documenta Ophthalmologica Proceedings Series, 49, 1987

[2] A. Heijl, G. Lindgren, J. Olsson, and P. Asman. \emph{Visual field interpretation with empiric probability maps}. Archives of Ophthalmology, 107, 1989
}

\seealso{ \code{\link{gloperc}} }

\author{Ivan Marin-Franch}