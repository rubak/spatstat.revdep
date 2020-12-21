\name{vfindexpmap}
\alias{vfindexpmap}
\title{probability values for visual field index}
\description{calculates the probability values for visual field index}
\usage{vfindexpmap( vfi )}
\arguments{
  \item{vfi}{ visual fields index }
}
\details{calculates the proability values for visual field index (see \code{\link{vfindex}})}

\references{
[1] B. Bengtsson and A. Heijl. \emph{A visual field index for calculation of glaucoma rate of progression}. American Journal of Ophthalmology, 145, 2008
}

\seealso{\code{\link{vfstats}}, \code{\link{vfstatspmap}}, \code{\link{vfindex}}}

\examples{
  vfi  <- vfindex( vf91016right )
  vfip <- vfindexpmap( vfi )
}

\author{Ivan Marin-Franch}