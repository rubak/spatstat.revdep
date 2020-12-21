\name{expectedPairs}
\alias{expectedPairs}
\alias{expectedCrossPairs}
\alias{expectedPairs_iso}
\alias{expectedCrossPairs_iso}
\alias{expectedPairs_kernloo}
\alias{expectedPairs_iso_kernloo}
\title{
Expected pairs in an inhomogeneous poisson process
}
\description{
Compute the expected number of pairs at a given displacement h in a poisson
process with a given intensity function. This corresponds to the integrals
\eqn{\gamma}{gamma} of Shaw et al. 2020. The various functions correspond
to the univariate and bivariate versions of the anisotropic or isotropic
versions of \eqn{\gamma}{gamma}. The final two options (\code{expectedPairs_kernloo}
and \code{expectedPairs_iso_kernloo}), provide implementations of the leave-out
kernel estimates of \eqn{\gamma}{gamma}: \eqn{\bar \gamma(h)}{gammabar(h)} and
\eqn{\bar \gamma ^\mathrm{iso}(r)}{gammabar^iso(r)}. In those cases, the point pattern
\code{X} itself is passed to the routine, rather than the (true or estimated) intensities
\code{rho} etc. The estimators for \eqn{\bar \gamma(h)}{gammabar(h)} are only applicable
to univariate processes. See Shaw et al, 2020 for details.}
\usage{
expectedPairs(rho, hx, hy=NULL, method=c("mc", "lattice"),
                tol=.005, dx=diff(as.owin(rho)$xrange)/200, maxeval=1e6,
                maxsamp=5e3)

expectedCrossPairs(rho1, rho2=NULL, hx, hy=NULL, method=c("mc", "lattice"),
                tol=.005, dx=diff(as.owin(rho1)$xrange)/200, maxeval=1e6,
                maxsamp=5e3)

expectedPairs_iso(rho, r, tol=.001, maxeval=1e6, maxsamp=5e3)

expectedCrossPairs_iso(rho1, rho2=NULL, r, tol=.001, maxeval=1e6, maxsamp=5e3)

expectedPairs_kernloo(X, hx,hy, sigma=bw.CvL, tol=.005, maxeval=1e6,
                            maxsamp=5e3, leaveoneout=TRUE)

expectedPairs_iso_kernloo(X, r, sigma=bw.CvL, tol=.001, maxeval=1e6,
                                    maxsamp=5e3, leaveoneout=TRUE)

}
\arguments{
  \item{rho1, rho2, rho}{
    Intensity functions, either of class \code{im} or \code{funxy}.  This may
    be produced by \code{\link{density.ppp}} or \code{\link{densityfun.ppp}},
    or provided by a fitted intensity model.
  }
  \item{X}{
    Point pattern of class \code{\link{ppp}} with the points of the pattern for which
    \eqn{\bar \gamma}{gammabar} is to be estimated.
  }
  \item{hx, hy}{
    For expectedPairs and expectedCrossPairs (i.e. \eqn{\gamma(h)}{gamma(h)}),
    the displacements \eqn{h \in \textrm{R}^2}{h in R^2} to evaluate \eqn{\gamma}{gamma}
    at. These can be in any format supported by xy.coords.
  }
  \item{r}{
    For the isotropic versions \eqn{\gamma^\mathrm{iso}(r)}{gamma[iso](r)}, the
    separations \eqn{r} at which \eqn{\gamma^\mathrm{iso}}{gamma[iso]} is to be
    evaluated.
  }
  \item{method}{
    Either mc (the default) or lattice. Compute integral using monte-carlo or
    on a lattice.
  }
  \item{tol}{
    A tolerance for how precise the integral should be. This is compared to a
    standard error for the mc estimate.
  }
  \item{sigma}{
    Smoothing bandwidth for direct kernel-based estimators \eqn{\bar \gamma}{gammabar}.
  }
  \item{leaveoneout}{
    Use leave-out estimators. This should generally be true except for the purpose of
    evaluating the bias of the standard estimators. See Shaw et al 2020 for details.
  }
  \item{maxeval}{
    Maximum number of evaluations of rho per iteration. Prevents memory-related
    crashes that can occur.
  }
  \item{maxsamp}{
    Maximum number of monte carlo samples per iteration. If this is too large,
    you may do more work than required to achieve tol.
  }
  \item{dx}{
    if method=="lattice", a lattice spacing for the computation. defaults to .01.
  }
}
\value{
  The return value is a numeric vector with length equal to the number of
  displacements h passed 
}
\references{
    T Shaw, J Møller, R Waagepetersen. 2020. \dQuote{Globally Intensity-Reweighted Estimators for
    \eqn{K}- and pair correlation functions}. arXiv:2004.00527 [stat.ME].
}
\author{
Thomas Shaw <shawtr@umich.edu>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pcfglobal}}, \code{\link{Kglobal}}, which use these functions to compute
the normalization functions \eqn{\gamma}{gamma}.
}
