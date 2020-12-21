\name{plot.locppm}
\alias{plot.locppm}
\alias{contour.locppm}
\title{
  Plot a Locally Fitted Poisson or Gibbs Model
}
\description{
  Plot an object of class \code{"locppm"}
  representing a locally-fitted Poisson or Gibbs point process model.
}
\usage{
\method{plot}{locppm}(x, \dots, what = "cg", which = NULL)

\method{contour}{locppm}(x, \dots, what = "cg", which = NULL)
}
\arguments{
  \item{x}{
    A locally-fitted Poisson or Gibbs point process model (object of class
    \code{"locppm"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link{plot.ssf}}
    to control the plot. 
  }
  \item{what}{
    What quantity to display. A character string.
    The default is to display the fitted coefficient vectors.
  }
  \item{which}{
    Which component(s) of the vector-valued quantity to display.
    An index or index vector.
  }
}
\details{
  These are methods for the generic commands \code{\link[graphics]{plot}}
  and \code{\link[graphics]{contour}}, for the class \code{"locppm"}.

  The argument \code{what} specifies what quantity will be displayed:
  \tabular{ll}{
    \code{"cg"} \tab Fitted coefficients of local model \cr
    \code{"vg"} \tab Local variance matrix for Gibbs model \cr
    \code{"vh"} \tab Local variance matrix for homogeneous model \cr
    \code{"tg"} \tab \eqn{t}-statistics based on \code{"coefs"} and \code{"vg"}
  }

  Typically these quantities are vector-valued (matrices are converted
  to vectors). The argument \code{which}, if present, specifies
  which elements of the vector are displayed. It may be any kind of
  index for a numeric vector.

  The plotting is performed by \code{\link{plot.ssf}}.
}
\value{
  \code{NULL}.
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
  \code{\link{locppm}}, 
  \code{\link{methods.locppm}}
}
\examples{
   fit <- locppm(swedishpines, ~1, sigma=9, nd=20,
             vcalc="hessian", locations="coarse")
   plot(fit)  
   plot(fit, what="Vg")  
}
\keyword{spatial}
\keyword{models}

