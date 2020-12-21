\name{methods.locmincon}
\alias{methods.locmincon}
\alias{as.ppp.locmincon}
\alias{print.locmincon}
\title{Methods for Local Cluster or Cox Models}
\description{
  Methods for various generic functions, for the class
  \code{"locmincon"} of locally fitted cluster or Cox point process models.
}
\usage{
  \method{as.ppp}{locmincon}(X, \dots)

  \method{print}{locmincon}(x, \dots)
}
\arguments{
  \item{x,X}{
    A locally-fitted Cox or cluster point process model (object of class
    \code{"locmincon"}). 
  }
  \item{\dots}{Additional arguments}
}
\value{
  \code{as.ppp} returns an object of class \code{"ppp"}.

  \code{print} returns \code{NULL}.
}
\details{
  Objects of class \code{"locmincon"} represent locally fitted
  cluster or Cox point process models.

  The functions documented here provided methods for this class,
  for the generic functions
  \code{\link[spatstat.geom]{as.ppp}}
  and
  \code{\link{print}}.
}
\references{
  \localpaper

  Baddeley, A., Rubak, E. and Turner, R. (2015)
  \emph{Spatial Point Patterns: Methodology and Applications with R}.
  Chapman and Hall/CRC Press.
}
\seealso{
  \code{\link{locmincon}}
}
\examples{
  example(locmincon)
  fit
  as.ppp(fit)  
}
\author{Adrian Baddeley}
\keyword{spatial}
\keyword{methods}

