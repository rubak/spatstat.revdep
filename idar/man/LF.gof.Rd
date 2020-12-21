\encoding{UTF-8}
\name{LF.gof}
\alias{LF.gof}

\title{Loosmore and Ford Goodness of Fit Test}
\description{ 
  Performs the Loosmore and Ford (2006) test or the Maximum Absolute Deviation test for a spatial point pattern.
}

\usage{
  LF.gof(X, rmin=NULL, rmax=NULL, na.rm=TRUE)
}

\arguments{
  \item{X}{An object resulting from the function \code{\link[spatstat]{envelope}},  i.e., with an attribute "\emph{simfuns}" (obtained using the argument  \code{savefuns=TRUE} in \code{\link[spatstat]{envelope}}) ,which is an object of class "\code{\link[spatstat]{fv}}" containing the summary functions computed for each of the simulated patterns.}
  \item{rmin}{Minimum value of the function argument r over which the maximum absolute deviation, or the integral, will be computed for the test.}
  \item{rmax}{Maximum value of the function argument r over which the maximum absolute deviation, or the integral, will be computed for the test.}
  \item{na.rm}{Should NA's be removed to compute the integral?}
 }

\details{
These function perform a tests for goodness-of-fit of a point pattern dataset to a point process model, based on Monte Carlo simulation from the model. The simulations should have been previously computed with the function \code{\link[spatstat]{envelope}}, applied with the argument \code{savefuns=TRUE} in order to save all the simulated functions, required for the computation of the test.

The test, popularized in the ecological field by Loosmore and Ford (2006) is also described in Diggle (2003, page 14), and according to Baddeley and Turner (2005) also in Diggle (1986) and  Cressie (1991, page 667, equation (8.5.42)).
If the arguments \code{rmin} and \code{rmax} are set to \code{NULL}, the integral of the GoF statistics will be computed over the complete range of r values.
}


\value{
  A list with the following components:
  \describe{
      \item{u}{The GoF statistic, i.e., the value of the integral over the range of \emph{r}'s}
      \item{p}{The p-value of the test}
      \item{na.count.by.r}{Number of \code{NA} values for each r. It helps to evaluate the reliability of the computed u's, specially for small \emph{r}'s}
   }
}

\seealso{
       \code{\link[spatstat]{dclf.test}} for an alternative implementation of the test in \pkg{spatstat}.
}


\references{
Cressie, N.A.C. (1991) \emph{Statistics for spatial data}. John Wiley and Sons, 1991. 

Diggle, P. J. (1986). Displaced amacrine cells in the retina of a rabbit : analysis of a bivariate spatial point pattern. \emph{J. Neuroscience Methods} 18, 115-125. 

Diggle, P.J. (2003) \emph{Statistical analysis of spatial point patterns}, Second edition. Arnold. 

Loosmore, N.B. and Ford, E.D. (2006) Statistical inference using the G or K point pattern spatial statistics. \emph{Ecology} 87, 1925-1931. 

}

\author{ Marcelino de la Cruz \email{marcelino.delacruz@urjc.es}}

\examples{
# some envelopesfor some idar function
# The argument "savefuns" must be set to "TRUE"
# BEWARE: in real tests nsim shpuld be 199 or higher
data(lansing)
maple.lambda<- density(unmark(lansing[lansing$marks=="maple"]))
multi.maple.sim <- multifocalsimulator(lansing, "maple", nsim=10,
                               simulate=expression(rpoispp(maple.lambda)))

isar.maple.env<- envelope(lansing, fun=isar, mimark="maple", nsim=10,
                                     simulate=multi.maple.sim, r=seq(0.01, 0.25, le=100),
                                      savefuns=TRUE)
# Estimate GoF test
LF.gof(isar.maple.env)

}

\keyword{ spatial }