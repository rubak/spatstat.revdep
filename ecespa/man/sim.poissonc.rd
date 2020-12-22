\name{sim.poissonc}
\alias{sim.poissonc}

\title{ Simulate Poisson Cluster Process }
\description{
  Generate a random point pattern, a simulated realisation of the Poisson Cluster Process
}
\usage{
sim.poissonc(x.ppp, rho, sigma)
}

\arguments{
  \item{x.ppp}{ Point pattern whose window and intensity will be simulated. An object with the
  \code{\link[spatstat.geom]{ppp}} format of \pkg{spatstat}. }
  \item{rho}{ Parameter \eqn{rho} of the Poisson Cluster process.  }
  \item{sigma}{ Parameter \eqn{sigma} of the Poisson Cluster process.  }
}
\details{
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

}
\value{The simulated point pattern (an object of class "\code{ppp}"). }
\references{ Diggle, P.J. 2003. \emph{Statistical analysis of spatial point patterns}. Arnold, London. }
\author{ Marcelino de la Cruz Rot }
\note{ This function can use the results of  \code{\link{pc.estK}} to simulate point patterns from a fitted model.
Be careful as the paramted returned by \code{\link{pc.estK}} is \eqn{sigma^2} while \code{sim.poissonc} takes 
its square root, i.e. \eqn{sigma}.
}
\section{Warning}{
This implementation simulates only point patterns within rectangular windows. Use  \code{\link{ipc.estK}} to fit and 
 \code{\link{rIPCP}} (or the \code{spatstat} functions) to simulate point patterns within irregular windows.
}
\seealso{ \code{\link{rIPCP}} to simulate inhomogeneous PCP; \code{\link[spatstat.core]{rNeymanScott}} 
and \code{\link[spatstat.core]{rThomas}} in \pkg{spatstat} }
\examples{


data(gypsophylous)

# set the number of simulations (nsim=199 or larger for real analyses)
nsim<- 39

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

## A kind of pointwise test of the pattern gypsophilous been a realisation
## of the fitted model, simulating with sim.poissonc and using function J (Jest).

gyps.env.sim <- envelope(gypsophylous, Jest,  nsim=nsim,
                    simulate=expression(sim.poissonc(gypsophylous,
		    sigma=sqrt(cosa.pc$sigma2), rho=cosa.pc$rho)))

plot(gyps.env.sim,  main="")


}
\keyword{spatial }

