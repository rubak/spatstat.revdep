\name{ttestmap}
\alias{ttestmap}
\title{
  Test of Effect in Locally Fitted Point Process Model
}
\description{
  Perform a local \eqn{t}-test for the presence of a covariate effect
  in a locally fitted Poisson or Gibbs point process model.
}
\usage{
ttestmap(object, term, \dots,
         method = c("exact", "hessian", "taylor"),
         grid = FALSE, 
         ngrid = NULL, grideps = NULL,
         verbose = TRUE)
}
\arguments{
  \item{object}{
    Locally fitted Poisson or Gibbs point process model
    (object of class \code{"locppm"}).
  }
  \item{term}{
    Term to be dropped from the model.
    A character string matching a term in the model formula
  }
  \item{\dots}{
    Ignored.
  }
  \item{method}{
    Choice of method to be used to evaluate the \eqn{t} statistic.
    See Details.
  }
  \item{grid}{
    Logical. If \code{FALSE}, compute the test statistic at all
    quadrature points of the model. If \code{TRUE}, compute 
    at a coarse grid of locations.
  }
  \item{ngrid}{
    Number of grid points (in each axis direction)
    for the coarse grid. Incompatible with \code{grideps}.
  }
  \item{grideps}{
    Spacing (horizontal and vertical) between grid points
    for the coarse grid. Incompatible with \code{ngrid}.
  }
  \item{verbose}{
    Logical value indicating whether to print progress reports.
  }
}
\details{
  The argument \code{object} should be a locally-fitted
  Poisson or Gibbs point process model (object of class
  \code{"locppm"} created by \code{\link{locppm}}).
  
  This function computes the local \eqn{t} test statistic
  for the test that a particular covariate effect in the model is zero.
  This is described in Baddeley (2016, sections 3 and 5).
}
\value{
  Object of class \code{"ssf"}.
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
 fit <- with(copper,
        locppm(Points, ~D, covariates=list(D=distfun(Lines)), nd=c(7,15)))
 plot(ttestmap(fit, "D"))
}
\keyword{spatial}
\keyword{htest}

