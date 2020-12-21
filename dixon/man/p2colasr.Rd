\name{p2colasr}
\alias{p2colasr}
\title{ P-value for a discrete distribution on small sample data }
\description{
   Computes the p-value  for a two-sided hypothesis test following
	Dixon's (2002:145) description of the method of Agresti & Min (2001).  
}
\usage{
p2colasr(Z, nsim = length(Z))
}
\arguments{
  \item{Z}{ \code{vector} with the observed Z-score in the first position and 
		all the simulated values behind. }
  \item{nsim}{ Number of simulated values. }
}
\value{P-value of the two-sided hypothesis test}
\references{ 
Agresti, A. & Min, Y. 2001. On small-sample confidence intervals 
for parameters in discrete distributions. \emph{Biometrics}, \bold{57}: 963-971.\cr\cr
Dixon, P.M. 2002. Nearest-neighbor contingency table analysis
 of spatial segregation for several species. \emph{Ecoscience}, \bold{9}(2): 142-151.\cr\cr
 }
\author{Marcelino de la Cruz Rot}
\note{ This function is usually not to be called by the user. It is internally used by \code{\link{dixon}}.}
\keyword{ distribution }
