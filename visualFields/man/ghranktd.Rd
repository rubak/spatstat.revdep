\name{ghranktd}
\alias{ghranktd}
\title{general height estimated from rank TD curve}
\description{
estimates the general height by comparison of the rank TD curve of the subject compared with mean normal rank TD curve
}
\usage{
ghranktd( td, minPts = 2, strategy = "isospaced", withinNormal = 95, pCentral = 1,
          link = make.link( "logit" ), scaleFactor = 52.4 )
}
\arguments{
  \item{td}{total-deviation values}
  \item{minPts}{minimum number of TD values within normal limits necessary to estimate the general height. Default is \code{2}. If there are less than 2 TD values within normal limits, then NA is returned}
  \item{strategy}{strategy for the comparison of the reconstructed within normal part of the rank TD curve with the mean normal rank TD curve. There are two posibilities at the moment \code{isospaced} and \code{parallel}. In the option \code{isospaced}, the TD values identified as within-normal limits are equally spaced and compared with the corresponding rank location of the mean normal TD curve. In the slightly more computationally complicated option \code{parallel}, the rank positions are obtained so that the all differences between subject's TD values and mean normal TD values in the corresponding rank position are equal. In both strategies, the mean rank position of the reconstructed rank TD curve corresponds with the mean rank position of the mean normal TD rank curve (for the pattern of locations 24-2 that is location 26.5 ). This way uniqueness is ensured on both strategies. Default is \code{parallel}. For more details see [1]}
  \item{withinNormal}{value or label used to identify the locations that are within normal limits. Default is \code{95}}
  \item{pCentral}{proportion of the TD rank curve to be used for the estimation of the general height. Default is \code{1}}
  \item{link}{link function to use with \code{\link{glm}} fit. See \code{\link{family}}. Defaulut is \code{logit}}
  \item{scaleFactor}{scale to use for ranks to make them be between 0 and 1. Default is \code{52.4}}
}

\value{returns the estimated general height for the test, or NA if the number of TD values within normal limits was less than \code{numPts}}

\seealso{\code{\link{pdval},\link{ghpostd}}}

\references{
[1] I. Marin-Franch, W. H. Swanson, and V. E. Malinovsky. A novel strategy for the
estimation of the general height of the visual field in patients with glaucoma.
Graefe's Archive Clinical Experimental Ophthalmology, 252(5):801-809, 2014.
}

\author{Ivan Marin-Franch}