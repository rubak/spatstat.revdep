\name{pcfglobal}
\alias{pcfglobal}
\alias{pcfcross.global}
\title{ (cross) pair correlation functions with a global intensity reweighting }
\description{
Compute \eqn{g_\textrm{global}}{g[global]} or \eqn{c_\textrm{global}}{c[global]}
}
\usage{
pcfglobal(X, lambda=NULL, ..., sigma=bw.CvL(X), r=NULL, rmax=NULL,
    kernel="epanechnikov", bw=NULL, stoyan=0.15, normtol=.005, ratio=FALSE,
    discrete.lambda=FALSE, divisor=c("r", "d"),
    leaveoneout=TRUE, interpolate=TRUE, interpolate.fac=10, exp_prs=NULL,
    interpolate.maxdx=diameter(as.owin(X))/100, dump=FALSE)

pcfcross.global(X,Y, lambdaX=NULL, lambdaY=NULL, ...,
    sigma=bw.CvL(X), r=NULL, rmax=NULL, kernel="epanechnikov", bw=NULL,
    stoyan=0.15, normtol=.005, ratio=FALSE, discrete.lambda=FALSE,
    divisor=c("r", "d"), analytical=NULL, interpolate=TRUE,
    interpolate.fac=10, exp_prs=NULL,
    interpolate.maxdx=diameter(as.owin(X))/100, dump=FALSE)

}
\arguments{
  \item{X, Y}{
    point process of type \code{\link{ppp}}, on which to evaluate the (cross) \eqn{K}-function
}
  \item{lambda, lambdaX, lambdaY}{
    intensity function estimates corresponding to \code{X} and \code{Y}. If omitted, intensity
    functions will be computed using \code{\link{density.ppp}} or \code{\link{densityfun.ppp}}
    (see discrete.lambda below)
}
  \item{\dots}{
    extra args passed to density.ppp or densityfun.ppp, if applicable.
}
  \item{sigma}{
    Bandwidth value to use for kernel-based intensity estimation, intensity functions and
    \code{exp_prs} are not provided by the user.
}
  \item{r}{
    Values of \eqn{r} to evaluate \eqn{K(r)} at. If omitted, a sensible default is chosen, using the same
    conventions as \code{\link{Kest}} and \code{\link{Kinhom}}.
}
  \item{rmax}{
    Maximum \eqn{r} to evaluate \eqn{K(r)} at. \code{rmax} is used to generate values for \code{r},
    if omitted. If missing, a sensible default is chosen.
}
  \item{kernel}{
    Kernel type for smoothing of pcf.
}
  \item{bw}{
    Kernel bandwidth for smoothing of pcf.
}
  \item{stoyan}{
    Coefficient for Stoyan's bandwidth selection rule. See \code{\link{pcf.ppp}}.
}
  \item{normtol}{
    A tolerance to use for expectedPairs or expectedCrossPairs when computing monte-carlo
    estimates of the normalizing factor \eqn{\gamma}{gamma}. Expressed as a maximum fractional
    standard error.
}
  \item{ratio}{
    If \code{TRUE}, assemble numerator and denominator of pcf estimator separately.
}
  \item{divisor}{
    Whether to use the evaluation distance (\code{"r"}) or the distance between points (\code{"d"})
    to normalize the contribution of each point pair.
}
  \item{analytical}{
    If \code{TRUE}, use Diggle-Jones weights %TODO: cite and comment about this.
}
  \item{discrete.lambda}{
    If \code{TRUE}, and intensity function(s) are not supplied, estimate intensities by
    interpolating the values on a discrete lattice (using \code{\link{interp.im}} and
    \code{\link{density.ppp}}), instead of exactly (using \code{\link{densityfun.ppp}}).
}
  \item{interpolate}{
    If \code{TRUE}, evaluate the expectedCrossPairs on a lattice and interpolate,
    rather than at the exact displacements observed in the pattern.
}
  \item{interpolate.fac}{
    If \code{interpolate}, the lattice spacing will be \code{sigma/interpolate.fac}.
}
  \item{leaveoneout}{
    Use the leave-one-out estimator for \eqn{\gamma}{gamma}. See Shaw et al 2020 for details.
}
  \item{exp_prs}{
    A function that returns values for
    \eqn{\gamma_\textrm{iso}(r)}{gamma[iso](r)}. If \eqn{\gamma}{gamma} is
    known explicitly, or the same calculation is being used for several
    point patterns, it can be much faster to compute it once and provide the function
    as \code{exp_prs}, since the computation of \eqn{\gamma}{gamma} is usually the
    slowest part.
}
  \item{interpolate.maxdx}{
    Upper bound on allowable lattice spacing for interpolation.
}
  \item{dump}{
    For debugging purposes, include computed values of \eqn{\gamma}{gamma} with the output,
    as \code{attr}s.
}
}
\value{
    The return value is an object of class \code{\link{fv}}, just as for \code{\link{pcf}}
    and \code{\link{pcfinhom}}. The object contains columns \code{r}, \code{theo}, and
    \code{global}, corresponding respectively to the argument \eqn{r}, the theoretical
    values of \eqn{g(r)} for a Poisson process, and \eqn{g_\mathrm{global}(r)}{g[global](r)}.
}
\references{
    T Shaw, J Møller, R Waagepetersen. 2020. \dQuote{Globally Intensity-Reweighted Estimators for
    \eqn{K}- and pair correlation functions}. arXiv:2004.00527 [stat.ME].
}
\author{
Thomas Shaw <shawtr@umich.edu>
}
\seealso{
\code{\link{expectedPairs}}
}
\examples{
rho <- funxy(function(x,y) 80*(1+x), owin())
X <- rpoispp(rho)
g <- pcfglobal(X)
#plot(g)
}
