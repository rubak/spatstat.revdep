\name{Smooth.locmincon}
\alias{Smooth.locmincon}
\title{
  Smooth a Locally Fitted Cluster or Cox Point Process Model
}
\description{
  Applies kernel smoothing to the fitted cluster parameters of a locally-fitted
  cluster or Cox point process model.
}
\usage{
\method{Smooth}{locmincon}(X, tau = NULL, ...)
}
\arguments{
  \item{X}{
    Object of class \code{"locmincon"}.
  }
  \item{tau}{
    Smoothing bandwidth.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link[spatstat]{Smooth.ppp}}
    controlling the smoothing and the pixel resolution.
  }
}
\details{
  An object of class \code{"locmincon"} represents
  a locally-fitted Cox or cluster point process model.
  It provides estimates of the cluster parameters at each of the
  data points of the original point pattern dataset.

  The parameter estimates will be smoothed using a Gaussian
  kernel with standard deviation \code{tau}.
}
\value{
  A pixel image or a list of pixel images.
}
\references{
  \localpaper

  Baddeley, A., Rubak, E. and Turner, R. (2015)
  \emph{Spatial Point Patterns: Methodology and Applications with R}.
  Chapman and Hall/CRC Press.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{locmincon}},
  \code{\link{Smooth.ppp}}
}
\examples{
   fit <- locmincon(redwood)
   Smooth(fit, tau=0.1)
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
