\name{pc.estK}
\alias{pc.estK}
\alias{Kclust}
\title{ Fit the Poisson Cluster Point Process by Minimum Contrast}
\description{
  Fits the Poisson Cluster point process to a point pattern dataset by the Method of Minimum Contrast. 
}
\usage{
pc.estK(Kobs, r, sigma2 = NULL, rho = NULL)
Kclust(r, sigma2, rho)
}
\arguments{
  \item{Kobs}{ Empirical \eqn{K}-function. }
  \item{r}{ Sequence of distances at which function \eqn{K} has been estimated. }
  \item{sigma2}{ Optional. Starting value for the parameter \eqn{sigma2} of the Poisson Cluster process. }
  \item{rho}{ Optional. Starting value for the parameter \eqn{rho} of the Poisson Cluster process. }
}
\details{
The algorithm fits the Poisson cluster point process to a point pattern, by finding the parameters of the Poisson cluster model
which give the closest match between the theoretical K function of the Poisson cluster process and the observed
K function. For a more detailed explanation of the Method of Minimum Contrast, see \code{\link[spatstat.core]{mincontrast}}
 in \pkg{spatstat} or Diggle (2003: 86). 

 The Poisson cluster processes are defined by the following postulates (Diggle 2003):
 \tabular{ll}{
         \emph{PCP1}\tab Parent events form a Poisson process with intensity \eqn{rho}.\cr
         \emph{PCP2}\tab Each parent produces a random number of offspring, according to a probability distribution \cr
	                    \tab \eqn{p[s]: s = 0, 1, 2, ...}\cr
         \emph{PCP3}\tab The positions of the offspring relative to their parents are distributed according to a bivariate pdf \eqn{h}.\cr
	 }
This implementation asumes that the probability distribution \eqn{p[s]} of offspring per parent is a Poisson distribution and 
that the position of each offspring relative to its parent follows a radially symetric Gaussian distribution with pdf

\deqn{h(x, y) = [1/(2*pi*sigma^2)]* exp[-(x^2+y^2)/(2*sigma^2)]}

The theoretical \eqn{K}-function of this Poisson cluster process is (Diggle, 2003):

\deqn{pi*r^2 + [1- exp(-r^2/4*sigma^2)]/rho}

The command \code{\link{Kclust}} computes the theoretical \eqn{K}-function of this Poisson cluster process and 
can be used to find some initial estimates of \eqn{rho} and \eqn{sigma^2}. In any case, the optimization usually finds the
correct parameters even without starting values for these parameters.

This Poisson cluster process can be simulated with \code{\link{sim.poissonc}}.

}
\note{The exponents \eqn{p} and \eqn{q} of the contrast criterion (see \code{\link[spatstat.core]{mincontrast}}) are fixed 
respectively to \eqn{p = 2} and \eqn{q = 1/4}. The \eqn{rmin} and \eqn{rmax} limits of integration of the 
contrast criterion are set up by the sequence of values of \eqn{r} and \eqn{Kobs} passed to \code{pc.estK}.}

\value{
    \item{sigma2}{Parameter \eqn{sigma^2}.}
    \item{rho }{Parameter \eqn{rho}. }
}
\references{ Diggle, P. J. 2003. \emph{Statistical analysis of spatial point patterns}. Arnold, London. }
\author{ Marcelino de la Cruz Rot, inspired by some code of  Philip M. Dixon }

\seealso{ \code{\link{ipc.estK}} for fitting the inhomogeneous Poisson cluster process; some functions in \pkg{spatstat}
( \code{\link[spatstat.core]{matclust.estK}} and \code{\link[spatstat.core]{lgcp.estK}}) fit other appropriate processes for clustered patterns;
\code{\link[spatstat.core]{mincontrast}} performs a more general implementation of the method of mimimum contrast.}
\examples{


data(gypsophylous)


# set the number of simulations (nsim=199 or larger for real analyses)
nsim<- 19

## Estimate K function ("Kobs").

gyps.env <- envelope(gypsophylous, Kest, correction="iso", nsim=nsim)

plot(gyps.env, sqrt(./pi)-r~r, legend=FALSE)

## Fit Poisson Cluster Process. The limits of integration 
## rmin and rmax are setup to 0 and 60, respectively. 

cosa.pc <- pc.estK(Kobs = gyps.env$obs[gyps.env$r<=60],
		           r = gyps.env$r[gyps.env$r<=60])

## Add fitted Kclust function to the plot.

lines(gyps.env$r,sqrt(Kclust(gyps.env$r, cosa.pc$sigma2,cosa.pc$rho)/pi)-gyps.env$r,
       lty=2, lwd=3, col="purple")

## A kind of pointwise test of the gypsophylous pattern been a realisation
## of the fitted model, simulating with sim.poissonc and using function J (Jest).

gyps.env.sim <- envelope(gypsophylous, Jest, nsim=nsim,
                    simulate=expression(sim.poissonc(gypsophylous,
		    sigma=sqrt(cosa.pc$sigma2), rho=cosa.pc$rho)))

plot(gyps.env.sim,  main="",legendpos="bottomleft")


}

\keyword{ spatial }

