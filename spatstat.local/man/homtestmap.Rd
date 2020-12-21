\name{homtestmap}
\alias{homteststat}
\alias{homtestmap}
\alias{update.homtestmap}
\title{
  Test Statistic for Homogeneity Test
}
\description{
  Compute the test statistic for the test of homogeneity
  of a locally-fitted Poisson or Gibbs point process model.
}
\usage{
homteststat(object, \dots, verbose = FALSE)

homtestmap(object, \dots,
           what=c("components", "statistic", "pvalue"),
           test = c("score", "taylor", "likelihood"),
           ladjust=c("none", "moment", "PSS"),
           calibrate=c("chisq", "Satterthwaite", "firstmoment"),
           simple = !is.null(theta0),
           theta0 = NULL,
           poolmoments=NULL,
           sigma = NULL, 
           saveall = FALSE, 
           use.fft = TRUE,
           verbose = TRUE)

\method{update}{homtestmap}(object, \dots, 
           what=NULL, test=NULL, ladjust=NULL,
           calibrate=NULL, saveall=FALSE, poolmoments=NULL)
}
\arguments{
  \item{object}{
    Locally-fitted point process (object of class \code{"locppm"})
    or an object previously computed by \code{homtestmap}.
  }
  \item{\dots}{
    For \code{homteststat}, arguments passed to \code{homtestmap}.
    For \code{homtestmap}, additional unmatched arguments are ignored.
  }
  \item{what}{
    Character string (partially matched)
    indicating whether to return the vector components
    of the local test statistic, or the value of the local test statistic, or
    the local \eqn{p}-values.
  }
  \item{test}{
    Character string (partially matched)
    indicating whether to perform
    the local score test (\code{test="score"}),
    or the local composite likelihood ratio test
    approximately (\code{test="taylor"}) or exactly 
    (\code{test="likelihood"}).
 }
  \item{ladjust}{
    Character string (partially matched)
    specifying an adjustment to the composite likelihood ratio test statistic.
  }
  \item{calibrate}{
    Character string (partially matched)
    specifying how to calculate \eqn{p}-values
    from the test statistic.
  }
  \item{simple}{
    Logical value indicating whether to treat the fitted model
    as a simple null hypothesis (\code{simple=TRUE}) or as an estimate
    of the parameters in a composite null hypothesis
    (\code{simple=FALSE}, the default).
  }
  \item{theta0}{
    Coefficient vector specifying a simple null hypothesis.
  }
  \item{poolmoments}{
    Logical value indicating how to calculate the reference distribution
    for the likelihood ratio test statistic (and thus how to
    calculate \eqn{p}-values). See Details.
  }
  \item{sigma}{
    Smoothing bandwidth.
  }
  \item{saveall}{
    Logical value indicating whether to compute a complete set of
    sufficient statistics and save them as an attribute of the result.
    See Details.
  }
  \item{use.fft}{
    For software testing purposes only.
    Logical value indicating whether to use
    data computed by the Fast Fourier Transform.
  }
  \item{verbose}{
    Logical value indicating whether to print progress reports.
  }
}
\details{
  These functions are used by \code{\link{homtest}} to
  perform a Monte Carlo test of the null hypothesis of
  homogeneity (i.e.\ constant parameter values) for the locally-fitted
  Poisson point process or Gibbs point process \code{object}.

  The function \code{homtestmap} computes
  either the local likelihood ratio test statistic
  or the local score test statistic.
  If \code{what="statistic"}, then the result is a scalar-valued
  function giving the local values of the test statistic.
  If \code{what="pvalue"}, the result is a scalar-valued function \eqn{p(v)}
  giving the local \eqn{p}-value at each location \eqn{v}.
  If \code{what="components"}, the result is a vector-valued
  function \eqn{T(v)} containing the components of the quadratic form;
  the squared norm of \eqn{T(v)} is
  equal to the desired test statistic at each location \eqn{v}.

  If \code{saveall=TRUE}, then a complete set of sufficient statistics is
  calculated and stored as an attribute of the result. This makes it
  possible to compute all of the statistics and \eqn{p} values
  described above. 

  The function \code{update.homtestmap}, a method for the
  generic function \code{\link[stats]{update}}, converts
  an object of class \code{"homtestmap"} from one of these formats to
  another, where possible. Except in trivial cases, this requires that
  the \code{"homtestmap"} object was computed with \code{saveall=TRUE}.
  
  The function \code{homteststat} computes the mean of
  the local test statistic or the mean
  of the local \eqn{p}-values over the
  observation window.

  To compute the \eqn{p}-values when \code{test="likelihood"}
  or \code{test="taylor"}, the values of the local likelihood ratio
  test statistic are referred to a gamma distribution whose first two
  moments are estimated from the data. If \code{poolmoments=FALSE},
  the local estimates of the moments are used; if
  \code{poolmoments=TRUE}, the spatial average of these estimates
  is used. The default is to use pooling whenever it is
  theoretically justified, namely when the template
  model is a stationary point process.
    
  Finer control over the computation is possible
  using the arguments \code{\dots} passed to \code{\link{locppm}}.
}
\value{
  For \code{homteststat}, a numeric value giving the test statistic.
  
  For \code{homtestmap} and \code{update.homtestmap},
  a spatially-sampled function object (class \code{"ssf"}; see
  \code{\link{ssf}}). This object also belongs to the special class
  \code{"homtestmap"} which has a print method.
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
  \code{\link{homtest}}
}
\examples{
   example(locppm)
   plot(H <- homtestmap(fit))
   H
}
\keyword{spatial}
\keyword{htest}
