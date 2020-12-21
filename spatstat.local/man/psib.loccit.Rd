\name{psib.loccit}
\alias{psib.loccit}
\alias{psib.locmincon}
\title{
  Sibling Probability of Locally Fitted Cluster Point Process
}
\description{
  Computes the sibling probability of a locally fitted
  cluster point process model.
}
\usage{
  \method{psib}{loccit}(object)

  \method{psib}{locmincon}(object)
}
\arguments{
  \item{object}{
    Fitted cluster point process model
    (object of class \code{"loccit"} or \code{"locmincon"}).
  }
}
\details{
  In a Poisson cluster process, two points are called \emph{siblings}
  if they belong to the same cluster, that is, if they had the same
  parent point. If two points of the process are
  separated by a distance \eqn{r}, the probability that
  they are siblings is \eqn{p(r) = 1 - 1/g(r)} where \eqn{g} is the
  pair correlation function of the process.
  
  The value \eqn{p(0) = 1 - 1/g(0)} is the probability that,
  if two points of the process are situated very close to each other,
  they came from the same cluster. This probability
  is an index of the strength of clustering, with high values
  suggesting strong clustering.

  This concept was proposed in Baddeley, Rubak and Turner (2015, page 479)
  and Baddeley (2016).

  The function \code{\link[spatstat.core]{psib}} is generic, with methods for
  \code{"kppm"}, \code{"loccit"} and \code{"locmincon"}.

  The functions described here are the methods for
  locally-fitted cluster models of class \code{"loccit"} and \code{"locmincon"}.
  They compute the spatially-varying sibling probability of the
  locally-fitted model.
}
\value{
  A spatially sampled function (object of class
  \code{"ssf"}) giving the spatially-varying sibling probability.
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
  \code{\link[spatstat.core]{psib}},
  \code{\link[spatstat.core]{kppm}},
  \code{\link{loccit}}
}
\examples{
   \dontrun{
   fit <- loccit(redwood, ~1, "Thomas")
   }
   \testonly{
   X <- redwood[owin(c(0,1), c(-1,-1/2))]
   fit <- loccit(X, ~1, "Thomas", nd=5, control=list(maxit=20))
   }
   fit
   plot(psib(fit))
}
\keyword{spatial}
\keyword{models}
