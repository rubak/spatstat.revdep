% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dm_functions.R
\name{controlmcmcSDA}
\alias{controlmcmcSDA}
\title{control.mcmcSDA}
\usage{
controlmcmcSDA(n.sim, burnin, thin, h, c1.h, c2.h)
}
\arguments{
\item{n.sim}{the number of iteration}

\item{burnin}{The number of burn-in}

\item{thin}{the number of thining}

\item{h}{tuning parameter of the proposal distribution used in the Langevin-Hastings MCMC algorithm (see Laplace.sampling); default is h=NULL and then set internally as 1.65/n(1/6), where n is the dimension of the random effect.}

\item{c1.h}{value of c1 used in the adaptive scheme for h; default is c1.h=0.01. See also 'Details' in PrevMap package}

\item{c2.h}{value of c2 used in the adaptive scheme for h; default is c2.h=1e-04. See also 'Details' in PrevMap package}
}
\value{
A list with processed arguments to be passed to the main function.
}
\description{
This function helps to define the number of iteration, burn-in, thining, and the tunning parameters of the adaptive MALA
}
\details{
To be used as one of the arguments of \code{\link{SDALGCPMCML}}
}
\examples{
n <- 545
h <- 1.65/(n^(1/6))
control.mcmc <- controlmcmcSDA(n.sim = 10000, burnin = 2000,
thin= 8, h=h, c1.h = 0.01, c2.h = 1e-04)
str(control.mcmc)
}
\seealso{
\link[PrevMap]{control.mcmc.MCML}
}
\author{
Olatunji O. Johnson \email{o.johnson@lancaster.ac.uk}

Emanuele Giorgi \email{e.giorgi@lancaster.ac.uk}

Peter J. Diggle \email{p.diggle@lancaster.ac.uk}
}
