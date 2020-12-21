\name{bw.locppm}
\alias{bw.locppm}
\title{
  Cross Validated Bandwidth Selection for Locally Fitted Point Process Model
}
\description{
  Uses cross-validation to select a smoothing bandwidth
  for locally fitting a Poisson or Gibbs point process model.
}
\usage{
bw.locppm(...,
          method = c("fft", "exact", "taylor"), 
          srange = NULL, ns = 9, sigma = NULL,
          additive = TRUE, 
          verbose = TRUE)
}
\arguments{
  \item{\dots}{
    Arguments passed to \code{\link[spatstat]{ppm}} to fit the homogeneous
    version of the model.
  }
  \item{method}{
    Method of calculation. The default \code{method="fft"} is much
    faster than the other choices. 
  }
  \item{srange}{
    Range of values of the smoothing parameter \code{sigma}
    to be searched. A numeric vector of length 2 giving the minimum
    and maximum values of \code{sigma}.
  }
  \item{ns}{
    Number of values of the smoothing parameter \code{sigma}
    in the range \code{srange} to be searched. A positive integer.
  }
  \item{sigma}{
    Vector of values of the smoothing parameter to be searched.
    Overrides the values of \code{ns} and \code{srange}.
  }
  \item{additive}{
    Logical value indicating whether to calculate the leverage
    approximation on the scale of the intensity (\code{additive=TRUE})
    or the log intensity (\code{additive=FALSE}).
    Applies only when \code{method = "taylor"}.
  }
  \item{verbose}{
    Logical value indicating whether to display progress reports.
  }
}
\details{
  This function determines the optimal value of the smoothing
  parameter \code{sigma} to be used in a call to \code{\link{locppm}}.

  The function \code{\link{locppm}} fits
  a Poisson or Gibbs point process model
  to point pattern data by local composite likelihood.
  The degree of local smoothing is controlled by a smoothing parameter
  \code{sigma} which is an argument to \code{\link{locppm}}.

  This function \code{bw.locppm} determines the optimal value of
  \code{sigma} by cross-validation.
  For each value of \code{sigma} in a search interval,
  the function \code{bw.locppm} fits the model locally
  with smoothing bandwidth \code{sigma}, 
  and evaluates the composite likelihood cross-validation criterion
  \code{LCV(sigma)} defined in Baddeley (2016), section 3.2.
  The value of \code{sigma} which maximises \code{LCV(sigma)} is returned.
}
\value{
  A numerical value giving the selected bandwidth.
  The result also belongs to the class \code{"bw.optim"}
  which can be plotted.
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
  Ns <- if(interactive()) 16 else 2
  b <- bw.locppm(swedishpines, ~1, srange=c(2.5,4.5), ns=Ns)
  b
  plot(b)
}
\keyword{spatial}
\keyword{models}
\keyword{smooth}
