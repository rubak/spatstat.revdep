\name{predict.locppm}
\alias{predict.locppm}
\alias{fitted.locppm}
\title{
  Prediction of a Locally Fitted Poisson or Gibbs Point Process Model
}
\description{
  Computes the fitted intensity of a locally-fitted Poisson point process
  model, or the fitted intensity, trend or conditional intensity
  of a locally-fitted Gibbs point process model.
}
\usage{
\method{fitted}{locppm}(object, \dots,
                        type = c("cif", "trend", "intensity"),
                        new.coef=NULL)

\method{predict}{locppm}(object, \dots,
                         type = c("cif", "trend", "intensity"),
                         locations=NULL, new.coef=NULL)
}
\arguments{
  \item{object}{
    A locally-fitted Poisson or Gibbs point process model (object of class
    \code{"locppm"}). 
  }
  \item{\dots}{
    Currently ignored.
  }
  \item{new.coef}{
    New vector or matrix of values for the
    model coefficients.
  }
  \item{locations}{
    Point pattern of locations where prediction should be computed.
  }
  \item{type}{
    Character string (partially matched) specifying the type of
    predicted value: the conditional intensity \code{"cif"} (the
    default), or the first order trend, or the intensity.
    For Poisson models all three options are equivalent.
  }
}
\details{
  These are methods for the generic functions
  \code{\link[stats]{fitted}}
  and 
  \code{\link[stats]{predict}}
  for the class \code{"locppm"} of  locally-fitted  Gibbs point process
  models.

  The \code{fitted} method computes,
  for each quadrature point \code{v}
  (or in general, at each point \code{v} where a local model was fitted),
  the intensity of the locally-fitted model at \code{v}.
  The result is a numeric vector.
  
  The \code{predict} computes the fitted intensity at any specified
  set of \code{locations}, and returns the result as an \code{ssf} object.
}
\value{
  For \code{fitted.locppm}, a numeric vector.

  For \code{predict.locppm}, an object of class \code{"ssf"}
  as described in \code{\link{ssf}}.
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
  \code{\link{locppm}}
}
\examples{
  fit <- locppm(cells, sigma=0.1, use.fft=TRUE)
  lam <- predict(fit)
}
\keyword{spatial}
\keyword{models}
