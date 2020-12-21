\name{predict.loccit}
\alias{predict.loccit}
\alias{fitted.loccit}
\title{
  Prediction for Locally-Fitted Cox or Cluster Model
}
\description{
  Computes the fitted intensity of a locally-fitted Cox process
  or cluster process model.
}
\usage{
  \method{predict}{loccit}(object, \dots)

  \method{fitted}{loccit}(object, \dots, new.coef=NULL)
}
\arguments{
  \item{object}{
    Locally fitted point process model (object of class \code{"loccit"}
    fitted by \code{\link{loccit}}).
  }
  \item{\dots}{
    Arguments passed to \code{\link{predict.locppm}}.
  }
  \item{new.coef}{
    New values for the fitted coefficients. A matrix
    in which each row gives the fitted coefficients at one of the
    quadrature points of the model.
  }
}
\details{
  The fitted intensity is computed.  
}
\value{
  An object of class \code{"ssf"} as described in \code{\link{ssf}}.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{loccit}}, \code{\link{predict.locppm}}.
}
\references{
  \localpaper

  Baddeley, A., Rubak, E. and Turner, R. (2015)
  \emph{Spatial Point Patterns: Methodology and Applications with R}.
  Chapman and Hall/CRC Press.
}
\examples{
  X <- redwood[owin(c(0,1), c(-1,-1/2))]
  fit <- loccit(X, ~1, "Thomas", nd=5, control=list(maxit=20))
  lam <- predict(fit)
}
\keyword{spatial}
\keyword{models}
\keyword{methods}
