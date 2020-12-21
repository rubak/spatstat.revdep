\name{vfstats}
\alias{vfstats}
\title{visual field stats}
\description{calculates the visual field stats}
\usage{vfstats( vf )}
\arguments{
  \item{vf}{\code{\link{vf}}-object with sensitivity thresholds}
}
\details{calculates the visual field stats: mean sensitivity (msens), std of sensitivities (ssens), mean total deviation (mtdev), std of total deviation (stdev), mean pattern deviation (mpdev), std of pattern deviation (spdev). All are weighted means and stds}

\references{
[1] A. Heijl, G. Lindgren, and J. Olsson. \emph{A package for the statistical analysis of visual fields}. Documenta Ophthalmologica Proceedings Series, 49, 1987
}

\seealso{\code{\link{vfstatspmap}}, \code{\link{vfindex}}, \code{\link{vfindexpmap}}}

\examples{
  vfs <- vfstats( vf91016right )
}

\author{Ivan Marin-Franch}