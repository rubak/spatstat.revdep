\name{vfstatspmap}
\alias{vfstatspmap}
\title{probability values for visual field global indices}
\description{calculates the probability values for visual field global indices}
\usage{vfstatspmap( vfindices )}
\arguments{
  \item{vfindices}{visual fields global indices}
}
\details{calculates the proability values for visual field global indices (see \code{\link{vfstats}})}

\references{
[1] A. Heijl, G. Lindgren, and J. Olsson. \emph{A package for the statistical analysis of visual fields}. Documenta Ophthalmologica Proceedings Series, 49, 1987
}

\seealso{\code{\link{vfstats}}, \code{\link{vfindex}}, \code{\link{vfindexpmap}}}

\examples{
  vfs  <- vfstats( vf91016right )
  vfsp <- vfstatspmap( vfs )
}

\author{Ivan Marin-Franch}