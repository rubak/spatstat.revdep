\name{plot.locmincon}
\alias{plot.locmincon}
\title{
  Plot a Locally Fitted Cluster or Cox Point Process Model
}
\description{
  Plot an object of class \code{"locmincon"}
  representing a locally-fitted cluster or Cox point process model.
}
\usage{
\method{plot}{locmincon}(x, \dots,
               how = c("exact", "smoothed"),
               which = NULL, sigma = NULL, do.points = TRUE)
}
\arguments{
  \item{x}{
    The model to be plotted.
    A locally-fitted cluster or Cox point process model (object of class
    \code{"locmincon"} or \code{"loccit"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat]{plot.ppp}}
    or \code{\link[spatstat]{plot.im}} to control the plot.
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
  \item{sigma}{
    Numeric. Smoothing bandwidth to be used if \code{how="smoothed"}.
  }
  \item{do.points}{
    Logical. Whether to display the original point data
    as well.
  }
} 
\details{
  This is a method for the generic command \code{\link[graphics]{plot}}
  for the class \code{"locmincon"}.

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
  \code{\link{locmincon}}, 
  \code{\link{methods.locmincon}}
}
\examples{
   X <- redwood[owin(c(0,1), c(-1,-1/2))]
   fitm <- locmincon(X, ~1, "Thomas", sigma=0.07)
   plot(fitm, how="smoothed")  
   plot(fitm, how="exact")  
}
\keyword{spatial}
\keyword{models}

