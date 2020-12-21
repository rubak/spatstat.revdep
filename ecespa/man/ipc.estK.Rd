\name{ipc.estK}
\alias{ipc.estK}
\alias{plot.ecespa.minconfit}
\alias{print.ecespa.minconfit}
\alias{ecespa.minconfit}
\title{ Fit the (In)homogeneous Poisson Cluster Point Process by Minimum Contrast }
\description{
  Fits the (In)homogeneous Poisson Cluster point process to a point pattern dataset by the Method of Minimum Contrast.
}
\usage{
ipc.estK(mippp, lambda = NULL, correction = "iso", r = NULL, sigma2 = NULL,
            rho = NULL, q = 1/4, p = 2)


\method{plot}{ecespa.minconfit}(x, type="L", add=FALSE, xlim=NULL, ylim=NULL, lwd=c(1,1),
           lty=c(1,2), col=c(1,2), main=NULL, ...)
}


\arguments{
  \item{mippp}{ Point pattern to which the  (I)PCP will be fitted. A point pattern with the \code{\link[spatstat.geom]{ppp}} format of \code{spatstat}. }
  \item{lambda}{ Optional. Values of the estimated intensity function as a pixel image (object of class "\code{\link[spatstat.geom]{im}}" of spatstat) giving the intensity values at all locations of \code{mippp}. }
  \item{correction}{ A character item selecting any of the options "border", "bord.modif", "isotropic", "Ripley" or
 "translate". It specifies the edge correction(s) to be applied in the computation of the \eqn{K(r)} function.  }
  \item{r}{  Numeric vector. The values of the argument \eqn{r} at which the \eqn{K(r)} functions  should be evaluated.  }
  \item{sigma2}{ Optional. Starting value for the parameter \eqn{sigma2} of the Poisson Cluster process. }
  \item{rho}{ Optional. Starting value for the parameter \eqn{rho} of the Poisson Cluster process. }
  \item{q}{ \eqn{q} exponent of the contrast criterion (see \code{\link[spatstat.core]{mincontrast}}). }
  \item{p}{ \eqn{p} exponent of the contrast criterion (see \code{\link[spatstat.core]{mincontrast}}).  }
  \item{x}{ An object of class 'ecespa.minconfit', resulting of  applying \code{ipc.estK} to fit a Poisson Cluster Process. }
  \item{type}{ Type of function to be ploted. If type="L", function \eqn{L(r)} [= \eqn{sqrt(K(r)/pi)-r}] is ploted. Otherwise, function \eqn{K(r)} is ploted. }
  \item{add}{ Logical. Should the curves be added to another plot? }
  \item{xlim}{ Vector setting the limits of the x-axis. }
  \item{ylim}{ Vector setting the limits of the y-axis. }
  \item{lwd}{ Vector (length=2) setting the line width for ploting the two functions. }
  \item{lty}{ Vector (length=2) setting the line type for ploting the two functions. }
  \item{col}{ Vector (length=2) setting the line color for ploting the two functions. }
  \item{main}{ Optional. Text to appear as a title of the plot. }
  \item{\dots}{ Additional graphical parameters passed to \code{link{plot}}. }
}
\details{
 The algorithm fits the (inhomogeneous) Poisson cluster point process (PCP) to a point pattern, by finding the parameters of the (inhomogeneous) Poisson cluster model
which give the closest match between the theoretical K function of the Poisson cluster process and the observed
K function. For a concise explanation of the PCP see \code{\link{pc.estK}}. For a more detailed explanation of the Method of Minimum Contrast, see \code{\link[spatstat.core]{mincontrast}}
 in \pkg{spatstat} or Diggle (2003: 86). 
 
 The inhomogeneous PCP can be thought of as a thinned process of an homogeneous PCP, where the spatially varying thinning probability
\eqn{f(s)} is related to the spatially varying intensity function \eqn{lambda(s)} as \eqn{f(s) = lambda(s)/max{lambda(s)}} (Waagepetersen, 2007).  
As the inhomogeneous K function  for  the IPCP coincides with the (homogeneous) K function for the corresponding homogeneous PCP, 
the parameters of the underlying homomgeneous PCP can be estimated as those that give the closest match between the theoretical K function 
for the homogeneous PCP and the empirical inhomogeneous K function  for  the observed IPCP.

This Poisson cluster process can be simulated with \code{\link{rIPCP}}.

}

\value{
\code{ipc.estK} gives an object of class '\code{ecespa.minconfit}', basically a list with the following components:
  
    		    
  \item{sigma2}{Parameter \eqn{sigma^2}.}
  \item{rho }{Parameter \eqn{rho}. }
  \item{d.theta }{Minimized value of the contrast criterion \eqn{D(theta)}. }
  \item{Kobs }{Values of the observed K-function. }
  \item{Kfit }{Values of the fitted K-function. }
  \item{r }{ Sequence of distances at which \code{Kobs} and \code{Kfit} have been estimated. }
  \item{data }{Original point pattern. }
  \item{lambda }{Original intensity function. }
  \item{dataname }{Name of the original point pattern. }
  \item{lambdaname }{Name of the original intensity function image. }
 \item{q}{ \eqn{q} exponent of the contrast criterion. }
  \item{p}{ \eqn{p} exponent of the contrast criterion.  }
 }
\references{
 Diggle, P. J. 2003. \emph{Statistical analysis of spatial point patterns}. Arnold, London.
 
Waagepetersen, R. P. 2007. An estimating function approach to inference for inhomogeneous Neymann-Scott processes. \emph{Biometrics} 63: 252-258. \url{https://doi.org/10.1111/j.1541-0420.2006.00667.x}.
 
 }
\author{ Marcelino de la Cruz Rot , inspired by some code of  Philip M. Dixon }

\seealso{some functions in \pkg{spatstat}: \code{\link[spatstat.core]{mincontrast}} for a more general implementation of the method of mimimum contrast;
\code{\link[spatstat.core]{matclust.estK}} and \code{\link[spatstat.core]{lgcp.estK}} fit other appropriate processes for clustered patterns.}

\examples{


#####################
## Same example as in pc.estK

data(gypsophylous)

# set the number of simulations (nsim=199 or larger for real analyses)
# and resolution for lambda (ngrid=200  or larger for real analyses)
nsim<- 19
ngrid <- 50

## Estimate K function ("Kobs").

gyps.env <- envelope(gypsophylous, Kest, correction="iso", nsim=nsim)

plot(gyps.env, sqrt(./pi)-r~r, legend=FALSE)

## Fit Poisson Cluster Process. The limits of integration 
## rmin and rmax are setup to 0 and 60, respectively. 

cosa.pc2 <- ipc.estK(gypsophylous, r = gyps.env$r[gyps.env$r<=60])

## Add fitted Kclust function to the plot.


plot(cosa.pc2, add=TRUE, lwd=c(3,3))


## A kind of pointwise test of the gypsophylous pattern been a realisation
## of the fitted model, simulating with rIPCP and using function J (Jest).


gyps.env.sim2 <- envelope(gypsophylous, Jest, nsim=nsim,
                    simulate=expression(rIPCP(cosa.pc2)))


plot(gyps.env.sim2,  main="",legendpos="bottomleft")

\dontrun{

#####################
## Inhomogeneous example

data(urkiola)
   
   #####################
   ## univariate case

   # get univariate pp
   I.ppp <- split.ppp(urkiola)$birch

   # estimate inhomogeneous intensity function
   I.lam <- predict (ppm(I.ppp, ~polynom(x,y,2)), type="trend", ngrid=ngrid)

   # Compute and plot envelopes to Kinhom, simulating from an Inhomogeneous
   #  Poisson Process:
   
   I2.env <- envelope( I.ppp,Kinhom, lambda=I.lam, correction="trans", 
                              nsim=nsim, simulate=expression(rpoispp(I.lam)))
   plot(I2.env, sqrt(./pi)-r~r, xlab="r (metres)", ylab= "L (r)", col=c(1,3,2,2),legend=FALSE) 

   # It seems that there is short scale clustering; let's fit an Inhomogeneous 
   # Poisson Cluster Process: 

   I.ki <- ipc.estK(mippp=I.ppp, lambda=I.lam, correction="trans")

   # Compute and plot envelopes to Kinhom, simulating from the fitted IPCP:

   Ipc.env <- Ki(I.ki, correction="trans", nsim=nsim, ngrid=ngrid)

   plot (Ipc.env, xlab="r (metres)", ylab= "L (r)")
 
   #####################
   ## bivariate case: test independence between birch and quercus in Urkiola

   J.ppp <- split.ppp(urkiola)$oak
   
   # We want to simulate oak from a homogeneous Poisson model:
   J.ppm <- ppm(J.ppp, trend=~1, interaction=Poisson() )
   
   IJ.env <- Kci (mod1=I.ki, mod2=J.ppm, nsim=nsim)
   
   plot(IJ.env, type=12)
   
   plot(IJ.env, type=21)
}


}

\keyword{ spatial }
