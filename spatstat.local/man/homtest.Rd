\name{homtest}
\alias{homtest}
\title{
  Homogeneity Test for Local Poisson or Gibbs Model
}
\description{
  Conducts a Monte Carlo test of homogeneity for a locally-fitted
  Poisson or Gibbs point process model.
}
\usage{
   homtest(X, \dots, nsim = 19,
           test = c("residuals", "score", "taylor", "likelihood"),
           locations = c("coarse", "fine", "split"),
           ladjust = NULL,
           use.fft = NULL,
           simul = NULL,
           verbose = TRUE, Xname = NULL)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
  }
  \item{\dots}{
    Additional arguments passed to \code{\link{locppm}} to determine the
    locally fitted model, and passed to \code{\link[spatstat]{ppm}} to
    determine the homogeneous model.
  }
  \item{nsim}{
    Number of simulations for the Monte Carlo test.
  }
  \item{test}{
    The local test statistic to be used:
    either \code{"likelihood"} for the local likelihood ratio test
    statistic, \code{"taylor"} for the Taylor approximation
    to the local likelihood ratio test statistic,
    \code{"score"} for the local score test statistic,
    or \code{"residuals"} for the squared local residuals.
  }
  \item{locations,use.fft}{
    Arguments passed to \code{\link{locppm}} to control the
    calculation of variances (if \code{method="local"}). See Details.
  }
  \item{ladjust}{
    Optional argument passed to \code{\link{homtestmap}} specifying a
    data-dependent adjustment of the test statistic.
  }
  \item{simul}{
    Optional information that determines 
    how to simulate the realisations from the null hypothesis.
    An expression in the \R language (that will be evaluated
    \code{nsim} times to obtain the simulated patterns), or a list
    that contains the simulated point patterns.
  }
  \item{verbose}{
    Logical value indicating whether to print progress reports.
  }
  \item{Xname}{
    Optional character string name for the dataset \code{X},
    to be printed in the test output.
  }
}
\details{
  This function performs a Monte Carlo test of the null hypothesis of
  homogeneity (i.e.\ constant parameter values) for the locally-fitted
  Poisson point process or Gibbs point process specified by the
  arguments.

  The type of test is controlled by the argument \code{test}.
  \itemize{
    \item
    \code{test="likelihood"}:
    the locally fitted model is computed as \code{locppm(X, \dots)}.
    The local composite likelihood ratio test statistic of this model
    is computed at each location,
    and the mean of this statistic over the window is computed.
    \item
    \code{test="taylor"}:
    the locally fitted model is computed as \code{locppm(X, \dots)}.
    The Taylor approximation to the
    local composite likelihood ratio test statistic of this model
    is computed at each location,
    and the mean of this statistic over the window is computed.
    \item
    \code{test="score"}:
    the locally fitted model is computed as \code{locppm(X, \dots)}.
    The local score test statistic of this model is computed at each location,
    and the mean of this statistic over the window is computed.
    \item
    \code{method="residuals"}:
    the homogeneous model is fitted as \code{ppm(X, \dots)}.
    The smoothed score residuals of this model are computed at each
    location, and the mean of the squared norm over the window
    is computed.
  }

  The test statistic is computed for the data pattern \code{X}
  and for each of \code{nsim} simulated realisations from the
  homogeneous model. The Monte Carlo \eqn{p}-value is computed.
}
\value{
  An object of class \code{"htest"} containing the test outcome.
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
  To compute the test statistic only, see
  \code{\link{homtestmap}}.
}
\examples{
   \dontrun{
   homtest(swedishpines)
   }
   \testonly{
   homtest(swedishpines, nd=5, nsim=5)
   }
}
\keyword{spatial}
\keyword{htest}

