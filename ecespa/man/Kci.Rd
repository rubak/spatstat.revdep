\name{Kci}
\alias{Kci}
\alias{Ki}
\alias{plot.ecespa.kci}
\alias{print.ecespa.kci}
\alias{ecespa.kci}

\title{ Test against non-Poisson (in-)homogeneous models }
\description{
  Functions to automate testing of 'univariate' and 'bivariate' point pattern hypothesis against non-Poisson (in-)homogeneous models.
}
\usage{
Kci(mod1, mod2, correction="trans", nsim=99, ngrid=200, nrep=1e+05,
     r=NULL, simu="both", spctype=1)

Ki(mod1, correction="trans", nsim=99, ngrid=200, nrep=1e+05, r=NULL,
    spctype=1)

## S3 ploth method for objects of class 'ecespa.kci':
\method{plot}{ecespa.kci}(x, type=1, q=0.025, kmean=TRUE, add=FALSE, maine=NULL, 
       xlabe=NULL, ylabe=NULL, xlime=NULL, ylime=NULL, 
       lty=c(1,2,3), col=c(1,2,3), lwd=c(1,1,1), ...)

}

\arguments{
  \item{mod1}{ Fitted model. An object of class \code{\link[spatstat.core]{ppm}} or  \code{\link{ecespa.minconfit}}. }
  \item{mod2}{ Fitted model. An object of class \code{\link[spatstat.core]{ppm}} or  \code{\link{ecespa.minconfit}}. }
  \item{correction}{ A character item selecting any of the options "border", "bord.modif",  or "translate". It specifies 
  the edge correction to be applied when computing K-functions. }
  \item{nsim}{ Number of simulated point patterns to be generated when computing the envelopes. }
  \item{ngrid}{ Dimensions (ngrid by ngrid) of a rectangular grid of locations where \code{\link[spatstat.core]{predict.ppm}} 
  would evaluate the spatial trend of the fitted models. }
  \item{nrep}{ Total number of steps (proposals) of Metropolis-Hastings algorithm that should be run by \code{\link[spatstat.core]{rmh}}
  to simulate models of class \code{ppm}.}
  \item{r}{ Numeric vector. The values of the argument \eqn{r} at which the \eqn{K(r)} functions  should be evaluated. }
  \item{simu}{ A character item indicating if both models will be simulated for the computation of the envelopes (\code{simu = "both"}) or 
  just the second model (\code{simu != "both"}). }
  \item{spctype}{ Type of 'pre-thinning' method employed by \code{\link{rIPCP}} in the simulation of
  \code{ecespa.minconfit} models. }
    
  \item{x}{ An object of class 'ecespa.kci'. The result of runing \code{Kci} or \code{Ki}. }
  \item{type}{ What to plot. One of 1 (\emph{K1}), 2 (\emph{K2}), 12 (\emph{K12}), 21 (\emph{K21}), 112 (\emph{K1-K12}) or 221 (\emph{K2-K21}). }
  \item{q}{Quantile for selecting the simulation envelopes.}
  \item{kmean}{Logical. Should the mean of the simulated envelopes be ploted?}
  \item{add}{Logical. Should the kci.plot be added to a previous plot? }
  \item{maine}{Title to add to the plot.}
  \item{xlabe}{Text  or expression to label the x-axis.}  
  \item{ylabe}{Text  or expression to label the y-axis. }
  \item{xlime}{Max and min coordinates for the x-axis.}
  \item{ylime}{Max and min coordinates for the y-axis.}
  \item{lty}{Vector with the line type for the estimated Kmm function, the simulated envelopes and the mean of the simulated envelopes. }
  \item{col}{Vector with the color for the estimated K-function, the simulated envelopes and the mean of the simulated envelopes.}
  \item{lwd}{Vector with the line width for the estimated K-function, the simulated envelopes and the mean of the simulated envelopes.}
  \item{\dots}{Additional graphical parameters passed to plot.}
}
\details{
  These functions are designed to automate the testing of 'univariate' and(/or) 'bivariate' point pattern hypotheses (based on K-functions) against non-Poisson (in-)homogeneous
  models. These non-Poisson (in-)homogeneous models should have been fitted with pseudolikelihood tools (\code{spatstat} \code{\link[spatstat.core]{ppm}} function)
or  with minimum contrast methods (\code{\link{ecespa.minconfit}}).

Function \code{Ki} is designed to test 'univariate' hypotheses. It will compute the (in-)homogeneous K-function (using \code{spatstat} \code{\link[spatstat.core]{Kinhom}} function) 
of the point pattern to which  the \code{\link[spatstat.core]{ppm}} or  \code{\link{ecespa.minconfit}} model has beeen fitted and will compute 'envelopes' simulating from the fitted model. 
The computed envelopes can be considered as a pointwise test of the point pattern been a realisation of the fitted model.

Function \code{Kci} is designed to test 'bivariate' hypotheses. It will compute the  (in-)homogeneous cross K-function 
(using \code{spatstat} \code{\link[spatstat.core]{Kcross.inhom}} function)  and will compute envelopes simulating from the fitted models.
As, when dealing with inhomogeneos patterns \eqn{K12} != \eqn{K21}, \code{Kci} will compute both functions. If \code{simu = "both"} (default option),
\code{Kci} will simulate from \code{mod2} to test K12 and from \code{mod1} to test K21. If \code{simu != "both"}, only \code{mod2} will be simulated. This option 
may be useful when only K12 is of interest. Function \code{Kci} will also compute univariate (in-) homogeneous K-functions and envelopes for each individual point pattern.

The S3 ploth method will plot the selected K-function and envelopes (actually, it will plot the most usual L-function = \eqn{sqrt[K(r)/pi]-r}). 
The appropriate K function can be selected with the argument \code{type}. If \code{type = 1} (default option), it will plot the univariate K function 
(of the analized model in \code{Ki} or of the first model [mod1] in \code{Kci}). If \code{type = 2}, it will plot the univariate K function of the second model 
(mod2 in \code{Kci}). When  \code{type = 12} or \code{type = 21}, it will plot respectively K12 or K21.
Options \code{type = 112} or \code{type = 221} will graph a kind of 'segregation test' (see \code{\link{K1K2}}), and will represent de differences 
K1-K12,  K2-K21 and their envelopes.

}
\value{
  Both \code{Kci} and \code{Ki} return an object of class \code{'ecespa.kci'}, basically a list with the following items:
  
  \item{r }{Numeric vector. The values of the argument \eqn{r} at which the \eqn{K(r)} functions  have been evaluated.}
  \item{kia }{Numeric vector. Observed (in-)homogeneous K function for \code{mod1} point pattern.}
  \item{kib }{ Numeric vector. Observed (in-)homogeneous K function for \code{mod2} point pattern.}
  \item{kci.ab.o }{ Numeric vector. Observed (in-) homogeneous cross K-function (K12) for \code{mod1} and \code{mod2} point patterns.}
  \item{kci.ba.o }{ Numeric vector. Observed (in-) homogeneous cross K-function (K21) for \code{mod2} and \code{mod1} point patterns.}
  \item{kci.ab.s }{ Matrix of simulated (in-) homogeneous cross K-function (K12) for \code{mod1} and \code{mod2} point patterns.}                
  \item{kci.ba.s }{ Matrix of simulated (in-) homogeneous cross K-function (K21) for \code{mod2} and \code{mod1} point patterns.}                
  \item{kib.s }{ Matrix of simulated (in-)homogeneous K function for \code{mod2} point pattern.}
  \item{kia.s }{ Matrix of simulated (in-)homogeneous K function for \code{mod1} point pattern.}
  \item{datanamea }{ Name of \code{mod1} point pattern.}
  \item{datanameb }{ Name of \code{mod2} point pattern.}
  \item{modnamea }{ Name of model \code{mod1}.}
  \item{modnameb }{ Name of model \code{mod2}.}
  \item{type }{ Type of analysis. "Kci" or "Ki".}
  
}
\references{

De la Cruz, M. and Escudero, A. 2008. Null models and tools for multivariate heterogeneous point patterns. \emph{Submitted}.

De Soto, L., Olano, J.M., Rozas, V. and De la Cruz, M. 2009. Release of \emph{Juniperus thurifera} woodlands from herbivore-mediated arrested succession in Spain. \emph{Applied Vegetation Science},  13: 1525. \url{http://dx.doi.org/10.1111/j.1654-109X.2009.01045.x}
}

\author{ Marcelino de la Cruz Rot }

\section{Warning }{As this implementation involves the use of images as the means of evaluation of the (inhomogeneous) spatial trend, and a mask based on
those images will be used as the point pattern window, the "Ripley's" or "isotropic" edge correction can not be employed.

It is usual that during the simulation process some warnings are produced. 
They are related to some simulated points being rejected as lying outside the specified window.

}

\note{ Even when one of the two point patterns is assumed to be homogeneous Poisson (and, apparently  not worth of fitting any model), an homogeneous
Poisson model can be easily fitted and passed to \code{Kci} with \code{\link[spatstat.core]{ppm}}. See the examples.
}

\examples{
  
    
    data(urkiola)
    
    # set the number of simulations (nsim=199 or larger for real analyses)
    nsim<- 5

  #####################
   ## univariate case

   # get univariate pp
   I.ppp <- split.ppp(urkiola)$birch

   # estimate inhomogeneous intensity function
   I.lam <- predict (ppm(I.ppp, ~polynom(x,y,2)), type="trend", ngrid=200)

   # Compute and plot envelopes to Kinhom, simulating from an Inhomogeneous
   #  Poisson Process:
   
   I2.env <- envelope( I.ppp,Kinhom, lambda=I.lam, correction="trans", 
                              nsim=nsim, simulate=expression(rpoispp(I.lam)))
   plot(I2.env, sqrt(./pi)-r~r, xlab="r (metres)", ylab= "L (r)", col=c(1,3,2,2),legend=FALSE) 

   # It seems that there is short scale clustering; let's fit an Inhomogeneous 
   # Poisson Cluster Process: 

   I.ki <- ipc.estK(mippp=I.ppp, lambda=I.lam, correction="trans")

   # Compute and plot envelopes to Kinhom, simulating from the fitted IPCP:

   Ipc.env <- Ki(I.ki, correction="trans", nsim=nsim, ngrid=200)

   plot (Ipc.env, xlab="r (metres)", ylab= "L (r)")
 
  \dontrun{
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


