\name{vfidefault}
\alias{vfidefault}
\docType{data}
\title{Settings of visual field index }
\description{Specifies the weighting to be applied at each location for averaging.}
\usage{data( vfidefault )}
\format{
vfisettings contains structures for different patterns of location testing. For instance, p24d2 or p30d2. Each structure contains 2 items:
  \describe{
    \item{\code{regweights}}{the actual weights to be applied to each region}
    \item{\code{locregions}}{mapping between locations in the visual-field testing pattern and regions with the same weight}
  }
}

\keyword{datasets}

\references{
[1] B. Bengtsson and A. Heijl. \emph{A visual field index for calculation of glaucoma rate of progression}. American Journal of Ophthalmology, 145, 2008
}

\author{Ivan Marin-Franch}