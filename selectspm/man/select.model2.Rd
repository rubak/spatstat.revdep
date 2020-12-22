\name{select.model2}
\alias{select.model2}
\alias{aic.function}
\alias{ipc.estK2}
\alias{envelope.selectedmod}
\alias{simulate.selectedmod}
\alias{plot.selectedmod}
\alias{print.selectedmod}
\title{
Fit and Select Point Pattern Models Based on Minimum Contrast and AIC
}
\description{
Fits Poisson, Poisson cluster, several inhomogeneous Poisson and several inhomogeneous Poisson cluster  processes to a spatial point pattern and select the best fitting based on AIC.
}
\usage{
select.model2(pp, sigmas, r, nlarge = 10000, q = 1/4, p = 2, correction = "trans")
ipc.estK2(mippp, lambda = NULL, correction = "iso", r = NULL, sigma2 = NULL, 
             rho = NULL, q = 1/4, p = 2, nlarge = NULL, ...)
aic.function(r, dtheta, npar)
  \method{envelope}{selectedmod}(Y,fun=NULL,nrank=1,nsim=99,dimyx=c(128,128),...)
  \method{simulate}{selectedmod}(object,nsim=99, seed=1,dimyx=c(128,128),...)
  \method{plot}{selectedmod}(x,...)
  \method{print}{selectedmod}(x,...)
  
}

\arguments{
  \item{pp}{
 Unmarked point pattern with the ppp format of spatstat
}
  \item{sigmas}{
Vector with the sigma values (standard deviations of the Gaussian kernel) for estimating different intensity surfaces with  density.ppp. 
}
  \item{r}{
Vector of values for the argument r at which the (in)homogeneous K function should be evaluated. First element should be 0.
}
  \item{nlarge}{
A number of points. In spatstat, if the number of points exceeds nlarge, then only the border correction will be computed (by default) for K(r). If you have a large number of points (n) and you want your border correction to be applied, set nlarge > n.
}
  \item{q}{
q parameter of the dtheta (i.e., minimum contrast) function.
}
  \item{p}{
p parameter of the dtheta (i.e., minimum contrast) function.
}
  \item{correction}{
Any selection of the options "border", "bord.modif", "isotropic", "Ripley", "translate", "translation", "none" or "best". It specifies the edge correction(s) to be applied.
}
\item{mippp}{
 Unmarked point pattern with the ppp format of spatstat
}
  \item{lambda}{
Optional. Values of the estimated intensity function. Either a vector giving the intensity values at the points of the pattern mippp, a pixel image (object of class "im") giving the intensity values at all locations, a fitted point process model (object of class "ppm" or "kppm") or a function(x,y) which can be evaluated to give the intensity value at any location.
}
  \item{sigma2}{
Starting  value in the optimization for the squared standard deviation of the Gaussian dispersion around parent points in the in(homogeneous) Poisson cluster process.
}
  \item{rho}{
Starting  value in the optimization for the intensity of parent points in the in(homogeneous) Poisson cluster process.
}
  \item{\dots}{
Additional arguments passed to the optim function (or to the plot, print, envelope and simulate methods).
}
\item{dtheta}{
Minimum contrast discrepancy, i.e., sum of squared diferences between the normalized empirical and theoretical K functions.
}
  \item{npar}{
Number of parameters fitted in the model.
}
 \item{x}{
An object of class "selectedmod", i.e., the result of using function select.model2().
}
\item{Y}{
An object of class "selectedmod", i.e., the result of using function select.model2().
}
\item{object}{
An object of class "selectedmod", i.e., the result of using function select.model2().
}
\item{fun}{
Function to compute the desired summary statistic for the point pattern.
}
\item{nsim}{
Number of simulated point patterns to be generated with simulate.selectedmod or when computing the envelopes with envelope.selectedmod.
}
\item{seed}{
A single value to set the random generator.
}

\item{nrank}{
Integer. Rank of the envelope value amongst the nsim simulated values. A rank of 1 means that the minimum and maximum simulated values will be used.
}
\item{dimyx}{
Pixel array dimensions to estimate the intensity of the point pattern when simulating from inhomogeneous models.
}

}
\details{
select.model2 is a wrap to fit and select different point processes using standard tools in of \pkg{spatstat} and of \pkg{ecespa}.
ipc.estK2 fits (in)homogeneous models as the function  \link[ecespa]{ipc.estK} of \pkg{ecespa} but, in addition, it allows seting the argument nlarge and passing options to \link[stats]{optim}.  AIC calculation (actually AICc) is made by aic.function. See more details in Jara et al. (2015).
}
\value{
select.model2 returns an object of class "\code{selectedmod}",i.e., a list with components:
\describe{
   \item{dtheta }{vector with the minimized discrepancy values for each fitted model.}
   \item{best.dtheta }{the minimum of the minimized discrepancy values.}
   \item{best.model }{ The best of the fitted models. }
   \item{models }{vector with the names of the fitted models}
   \item{HPPs }{A list with the intensity objects employed to fit each inhomogeneous model}
   \item{sigmas }{vector with the sigma values (standard deviations of the Gaussian kernel) employed to estimate the  different intensity surfaces. }
   \item{aics }{vector with the AIC values of each model.}
   \item{Kas }{data.frame with the empirical K functios employed to fit each model.}
   \item{pp }{ analyzed point pattern, with the ppp format of spatstat.}
}	  

\code{ipc.estK} gives an object of class '\code{ecespa.minconfit}', basically a list with the following components:
\describe{ 
 \item{sigma2}{Parameter \eqn{sigma^2}.}
  \item{rho }{Parameter \eqn{rho}. }
  \item{d.theta }{Minimized value of the contrast criterion \eqn{D(theta)}. }
  \item{Kobs }{Values of the observed K-function. }
  \item{Kfit }{Values of the fitted K-function. }
  \item{r }{ Sequence of distances at which \code{Kobs} and \code{Kfit} have been estimated. }
  \item{data }{Original point pattern. }
  \item{lambda }{Original intensity object. }
  \item{dataname }{Name of the original point pattern. }
  \item{lambdaname }{Name of the original intensity object. }
 \item{q}{ \eqn{q} exponent of the contrast criterion. }
  \item{p}{ \eqn{p} exponent of the contrast criterion.  }
}
\code{aic.function} returns a one-row data.frame, with the following items:
  \describe{
  \item{n}{number of observations employed to compute AIC; i.e., number of r values where K(r) was estimated.}
  \item{K}{number of parameters of the model.}
  \item{RSS}{"Resiidual Sum of Squares". It is the value of the discrepacy function Dtheta.}
  \item{LL}{Loglikelihood.}
  \item{AIC}{AIC value.}
  \item{AICc}{Small sample AIC value.}
}
\code{envelope.selectedmod} returns a an object of class "fv", see \link[spatstat.core]{fv.object}, which can be printed and plotted directly. Essentially a data frame containing columns 
  \describe{
  \item{r}{the vector of values of the argument r at which the summary function fun has been estimated.}
  \item{obs}{values of the summary function for the data point pattern.}
  \item{lo}{lower envelope of simulations.}
  \item{hi}{upper envelope of simulations.}
  \item{mmean}{estimated theoretical value of the summary function, computed by averaging simulated values.}
 }
 \code{simulate.selectedmod} returns a list of \emph{nsim} point patterns (with the \link[spatstat.geom]{ppp} format of spatstat, simulated according to the best fitted model.
 
  
}
\references{
Jara, A., De la Cruz, M., Espinosa, C.I., Mendez, M. & Escudero, A. (2015). Does spatial heterogeneity blur the signature of dispersal syndromes on spatial patterns of woody species? A test in a tropical dry forest. \emph{Oikos}. http://dx.doi.org/10.1111/oik.02098.
}
\author{
Marcelino de la Cruz
}

\examples{
\dontrun{
# Get the data
data(lansing)

# Split the multivariate pp in their individual components
lansing.sp<-split(lansing)

# Define the sequence of r's at which estimate K(r)
r<-  seq(0,0.25,le=101)

# Define different standard deviations for the Gaussian kernel 
# to estimate different intensity surfaces
sigmas<- seq(0.1,1,by=0.05)

# Note that lansing is defined in a (0,1) x (0,1) window and this affects
# the election of r and sigma values

# Fit 40 models (1 Poisson, 1 Poisson cluster, 19 inhomogeneous Poisson
#  and 19 inhomogeneous Poisson cluster)  to maple and select the better ones
maple.model <- select.model2(lansing.sp$maple, sigmas=sigmas, r=r)
# show the AICc value and the fitted parameters for the best model in each class
maple.model

# Draw the empirical and theoretical models to visually asses the fitting.
# P = Poisson; HPP= heterogeneous (i.e. inhomogeneous) Poisson;
# PC = Poisson cluster; HPC=heterogeneous (i.e. inhomogeneous) Poisson cluster
plot(maple.model)

# Compute and plot envelopes for the K function according to the best fitted model.
maple.env <- envelope(maple.model, nsim=19)
plot(maple.env, sqrt(./pi)-r~r, legend=F)

# simulate 10 point patterns according to the best fitted model

maple.simu <- simulate(maple.model, nsim=10)

maple.simu

# FIt and select models to all species
lansing.models<-lapply(lansing.sp, function(x) select.model2(x, sigmas=sigmas, r=r))
lapply(lansing.models, function(x) x)
}
}

\keyword{ spatial }
