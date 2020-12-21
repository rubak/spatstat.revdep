\name{rIPCP}
\alias{rIPCP}
\title{ Simulate Inhomogeneous Poisson Cluster Process }
\description{
  Generate a random point pattern, a simulated realisation of the Inhomogeneous Poisson Cluster Process. 
}
\usage{
rIPCP(x, lambda = NULL, type = 1, lmax = NULL, win = owin(c(0, 1), c(0, 1)), ...)
}
\arguments{
  \item{x}{ an object of class \code{'ecespa.minconfit'}, resulting from the function \code{\link{ipc.estK}}. }
  \item{lambda}{ Optional. Values of the estimated intensity function as a pixel image (object of class "im" of \code{spatstat}) giving the intensity values at all locations. }
  \item{type}{ Type of 'prethining' employed in the simulation. See details. }
  \item{lmax}{ Optional. Upper bound on the values of lambda. }
  \item{win}{ Optional. Window of the simulated pattern. }
  \item{\dots}{ Optional. Arguments passed to \code{\link[spatstat.geom]{as.im}}. }
}
\details{
  This function simulates the Inhomogeneous Poisson Cluster process from an object of class \code{'ecespa.minconfit'}, resulting from fitting an IPCP to some 'original' point pattern
  using the function \code{\link{ipc.estK}}.   Following the approach of Waagepetersen (2007), the simulation involves a first step in which an homogeneous aggregated pattern 
  is simulated (from the fitted parameters  of the \code{'ecespa.minconfit'} object, using function \code{\link[spatstat.core]{rThomas}}   of \code{spatstat}) and a second one in which
  the homogeneous pattern is thinned with a spatially varying thinning probability \emph{f (s)} proportional to the spatially varying intensity, i.e. \emph{f (s) = lambda(s) / max[lambda(s)]}.
  To obtain a 'final' density similar to that of the original point pattern, a "prethinning" must be performed. There are two alternatives. If the argument \code{'type'} is set equal to '1', 
  the expected number of points per cluster (\emph{mu} parameter of \code{\link[spatstat.core]{rThomas}} is thinned as \emph{mu <- mu.0 / mean[f(s)]}, where \emph{mu.0} is the 
  mean number of points per cluster of the original pattern. This alternative produces point patterns most similar to the 'original'.  If the argument \code{'type'} is set equal to '2', 
  the fitted intensity of the Poisson process of cluster centres (\emph{kappa} parameter of \code{\link[spatstat.core]{rThomas}}, i.e. the intensity of 'parent' points) is thinned 
  as \emph{kappa <- kappa / mean[f(s)]}. This alternative produces patterns more uniform than the 'original' and it is provided only for experimental purposes.
}
\value{
  A point pattern, with the format of the \code{\link[spatstat.geom]{ppp}} objects of \code{spatstat}.
}
\references{ 
Waagepetersen, R. P. 2007. An estimating function approach to inference for inhomogeneous Neymann-Scott processes. \emph{Biometrics} 63: 252-258. \url{https://doi.org/10.1111/j.1541-0420.2006.00667.x}.

}
\author{ Marcelino de la Cruz Rot }

\seealso{ \code{\link{sim.poissonc}} to simulate homogeneous PCP; \code{\link[spatstat.core]{rNeymanScott}} 
and \code{\link[spatstat.core]{rThomas}} in \pkg{spatstat} are the basis of this function }
\examples{

  
    data(gypsophylous)
  
    plot(gypsophylous) 
    
    ## It 'seems' that the pattern is clustered, so 
    ## fit a Poisson Cluster Process. The limits of integration 
    ## rmin and rmax are setup to 0 and 60, respectively.
    
   cosa.pc2 <- ipc.estK(gypsophylous, r = seq(0, 60, by=0.2))

   ## Create one instance of the fitted PCP:

   pointp <- rIPCP( cosa.pc2)
   
   plot(pointp)
   
   
   
    #####################
    ## Inhomogeneous example

    data(urkiola)

    # get univariate pp
    I.ppp <- split.ppp(urkiola)$birch

    plot(I.ppp)

    #estimate inhomogeneous intensity function
    I.lam <- predict (ppm(I.ppp, ~polynom(x,y,2)), type="trend", ngrid=200)

    # It seems that there is short scale clustering; lets fit an IPCP: 

    I.ki <- ipc.estK(mippp=I.ppp, lambda=I.lam, correction="trans")

    ## Create one instance of the fitted PCP:

    pointpi <- rIPCP( I.ki)
   
    plot(pointpi)
 

}
\keyword{ spatial }
