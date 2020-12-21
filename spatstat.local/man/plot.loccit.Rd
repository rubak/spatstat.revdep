\name{plot.loccit}
\alias{plot.loccit}
\title{
  Plot a Locally Fitted Cluster or Cox Point Process Model
}
\description{
  Plot an object of class \code{"loccit"}
  representing a locally-fitted cluster or Cox point process model.
}
\usage{
\method{plot}{loccit}(x, \dots,
               what = c("modelpar", "coefs", "lambda"),
               how = c("smoothed", "exact"), which = NULL,
               pre=NULL, post=NULL)
}
\arguments{
  \item{x}{
    The model to be plotted.
    A locally-fitted cluster or Cox point process model (object of class
    \code{"loccit"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat]{plot.ppp}}
    or \code{\link[spatstat]{plot.im}} to control the plot.
  }
  \item{what}{
    Character string determining which quantities to display:
    \code{"modelpar"} for the cluster model parameters,
    \code{"coefs"} for the trend coefficients,
    or \code{"lambda"} for the fitted intensity.
  }
  \item{how}{
    Character string determining whether to display the
    fitted parameter values at the data points (\code{how="exact"})
    or the smoothed fitted parameters as pixel images (\code{how="smoothed"}).
  }
  \item{which}{
    Optional. Which component(s) of the vector-valued quantity to display.
    An index or index vector. Default is to plot all components.
  }
  \item{pre,post}{
    Transformations to apply before and after smoothing.
  }
} 
\details{
  This is a method for the generic command \code{\link[graphics]{plot}}
  for the class \code{"loccit"}.

  The argument \code{which}, if present, specifies
  which fitted parameters are displayed. It may be any kind of
  index for a numeric vector.

  The quantities are computed at irregularly-placed points.
  If \code{how="exact"} the exact computed values
  will be displayed as circles centred at the locations where they
  were computed. If \code{how="smoothed"} these
  values will be kernel-smoothed using \code{\link[spatstat]{Smooth.ppp}}
  and displayed as a pixel image.  
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
  \code{\link{loccit}}
}
\examples{
   X <- redwood[owin(c(0,1), c(-1,-1/2))]
   fitc <- loccit(X, ~1, "Thomas", nd=5, control=list(maxit=20))
   plot(fitc, how="exact")  
   plot(fitc, how="smoothed")  
}
\keyword{spatial}
\keyword{models}

