\name{rrmix}
\alias{rrmix}
\alias{rrim}
\title{Spatial relative risk surface generation}
\description{
Generates an appropriately scaled spatial (bivariate) relative risk surface using a supplied control density and \eqn{N} isotropic Gaussian-style hotspots.
}
\usage{
rrmix(g, rhotspots, rsds, rweights, rbase = 1, log = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{
  A pixel \code{\link[spatstat.geom]{im}age} representing the control density; this will be internally rescaled to integrate to 1 if it does not already do so.
  }
  \item{rhotspots}{
  A \eqn{2 \times N} matrix giving the centers of the \eqn{N} peaks and troughs in the relative risk density.
  }
  \item{rsds}{
  A positive numeric vector of length \eqn{N} giving the isotropic standard deviations for each relative Gaussian peak or trough.
  }
    \item{rweights}{
  A vector of length \eqn{N} giving relative weightings for each peak (positive weights) or trough (negative).
  }
%    \item{window}{
%  An object of class \code{\link[spatstat.geom]{owin}} giving the observational window on which the surfaces are defined. If \code{NULL} (default), the polygonal version of \code{Window(g)} is used.
%  }
  \item{rbase}{
  The base level of the relative risk surface (default is 1). The peaks and troughs will be added or subtracted from this base level prior to normalisation.
  }

  \item{log}{
  A logical value. If \code{TRUE} (default), the relative risk surface is returned logged.
  }
}
\details{
A useful tool for the comparison of two estimated density functions on the same spatial region \eqn{W \subset R^2} is the \emph{relative risk function}, \eqn{r}, (Bithell, 1990; 1991; Kelsall and Diggle, 1995), defined simply as a density-ratio:
\deqn{
r(x)=f(x) / g(x);   x \in W.
}
Various methods have been developed to improve estimation of \eqn{r}, most commonly with a motivation in geographical epidemiology, where the `numerator' density \eqn{f} pertains to the observed disease cases and the `denominator' density \eqn{g} reflects the distribution of the at-risk controls (Kelsall and Diggle, 1995; Hazelton and Davies, 2009; Davies and Hazelton, 2010). To test newly developed methodology, simulations based on known relative risk scenarios are usually necessary. This function allows the user to design such scenarios, as used in Hazelton and Davies (2009), Davies and Hazelton (2010), and Davies (2013) for example.

This function calculates a relative risk surface based on \eqn{N} Gaussian-style `bumps' added and subtracted from a base level of \code{rbase}, with the peaks and troughs centered at the coordinates given by \code{rhotspots} with relative weights of \code{rweights} and isotropic standard deviations of \code{rsds}. The risk surface \eqn{r} is computed as

\eqn{r(x) \propto } \code{rbase} \eqn{+ \sum_{i=1}^{N}} \code{rweights[}\eqn{i}\code{]}\eqn{*exp(-0.5*}\code{rsds[}\eqn{i}\code{]}\eqn{^(-2)*}||\eqn{x-}\code{rhotspots[,}\eqn{i}\code{]}||\eqn{^2)}

where || . || denotes Euclidean norm. Because \eqn{f} and \eqn{g} are both densities, the risk surface as defined above must then be rescaled with respect to the supplied control density \eqn{g} (argument \code{g}) to ensure that
\deqn{
\int_W r(x)g(x) dx = 1
}
This is automatically performed inside the function. The case density that gives rise to the designed \eqn{r} is then easily recovered because \eqn{f = r * g}. By default, the function returns the log-relative risk surface \eqn{\log r = \log f - \log g} alongside the case and control densities.
}

\value{
An object of class \code{rrim}. This is a \code{\link[spatstat.geom]{solist}} of three pixel \code{\link[spatstat.geom]{im}}ages: \code{f} as the case density, \code{g} the control density (a copy of the argument of the same name, integrating to 1), and \code{r} as the (log) relative risk surface.
}

\references{
Bithell, J.F. (1990), An application of density estimation to geographical epidemiology, \emph{Statistics in Medicine}, \bold{9}, 691-701.\cr\cr
Bithell, J.F. (1991), Estimation of relative risk functions, \emph{Statistics in Medicine}, \bold{10}, 1745-1751.\cr\cr
Davies, T.M. (2013), Jointly optimal bandwidth selection
for the planar kernel-smoothed density-ratio, \emph{Spatial and
Spatio-temporal Epidemiology}, \bold{5}, 51-65.\cr\cr
Davies, T.M. and Hazelton, M.L. (2010), Adaptive kernel estimation of spatial relative risk, \emph{Statistics in Medicine}, \bold{29}(23) 2423-2437.\cr\cr
Kelsall, J.E. and Diggle, P.J. (1995a), Kernel estimation of relative risk, \emph{Bernoulli}, \bold{1}, 3-16.\cr\cr
}

\author{A.K. Redmond and T.M. Davies}

\examples{
set.seed(1)
gg <- rgmix(3,window=toywin,S=matrix(c(0.08^2,0,0,0.1^2),nrow=2),p0=0.2)

rho <- rrmix(g=gg,
             rhotspots=cbind(c(0.8,0.3),c(0.4,0.4),c(0.6,0.5),c(0.3,0.5)),
             rsds=c(0.005,0.025,0.01,0.025),
             rweights=c(3,2,10,5)*10)

rho.sample <- rrpoint(c(400,800),rho,toywin)


par(mfrow=c(2,2))
plot(rho$g,main="control density")
plot(rho$f,main="case density")
plot(rho$r,main="log relative risk surface")
plot(rho.sample$controls,main="sample data")
points(rho.sample$cases,col=2)
legend("topright",col=2:1,legend=c("cases","controls"),pch=1)
}


