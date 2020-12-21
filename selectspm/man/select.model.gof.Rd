\name{select.model.gof}
\alias{select.model.gof}
\alias{envelope.selectedmodgof}
\alias{simulate.selectedmodgof}
\alias{plot.selectedmodgof}
\alias{print.selectedmodgof}
\title{
Fit and Select Point Pattern Models Based on Minimum Contrast and Goodness of Fit
}
\description{
Fits Poisson, Poisson cluster, several inhomogeneous Poisson and several inhomogeneous Poisson cluster  processes to a spatial point pattern and select the best fitting based on goodness of fit.
}
\usage{
select.model.gof(pp, sigmas, r, nlarge = 10000, q = 1/4, p = 2, correction = "trans",
                          sigma2=NULL, rho=NULL, lower=NULL, upper=NULL, parscale=c(1,1),
	                   dimyx=c(128,128), nsim=99, seed=1, correct.lambda=10) 
  \method{plot}{selectedmodgof}(x,...)
  \method{print}{selectedmodgof}(x,...)
  \method{envelope}{selectedmodgof}(Y,fun=NULL,nrank=1,nsim=99,dimyx=c(128,128),...)
  \method{simulate}{selectedmodgof}(object, nsim=99, seed=1, dimyx=c(128,128),...)
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
\item{sigma2}{
Starting  value in the optimization for the squared standard deviation of the Gaussian dispersion around parent points in the in(homogeneous) Poisson cluster process.
}
  \item{rho}{
Starting  value in the optimization for the intensity of parent points in the in(homogeneous) Poisson cluster process.
}
 \item{lower}{
Vector of length two with the lowest allowed values for sigma2 and rho.
}
\item{upper}{
Vector of length two with the largest allowed values for sigma2 and rho.
}
\item{parscale}{
Initial values for sigma2 and rho for one of  the optimization approaches .
}
\item{dimyx}{
Pixel array dimensions to estimate the intensity of the point pattern when simulating from inhomogeneous models.
}
\item{nsim}{
Number of simulated point patterns to be generated for the estimation of the GoF of each fitted model, or when using  simulate.selectedmodgof or  envelope.selectedmodgof.
}
\item{seed}{
A single value to set the random generator.
}
\item{correct.lambda}{
Fraction of the lowest positive  intensity value that will be employed to replace intensity estimates == 0. 
}
  \item{\dots}{
Additional arguments passed to  the plot, print, envelope and simulate methods.
}
 \item{x}{
An object of class "selectedmodgof", i.e., the result of using function select.model.gof().
}
 \item{Y}{
An object of class "selectedmodgof", i.e., the result of using function select.model.gof().
}
\item{object}{
An object of class "selectedmodgof", i.e., the result of using function select.model.gof().
}
\item{fun}{
Function to compute the desired summary statistic for the point pattern.
}
\item{nrank}{
Integer. Rank of the envelope value amongst the nsim simulated values. A rank of 1 means that the minimum and maximum simulated values will be used.
}
}

\details{
select.model.gof is a wrap to fit and select different point processes using standard tools in \pkg{spatstat} and in \pkg{ecespa}. The basic framework consists in choosing among the provided \emph{sigma} values the bandwith which produces the best fitting of an inhomogeneous Poisson and the one which produces the best fitting of an inhomogeneous Poisson cluster process. The goodness of fit (GoF) of these models are  compared to the GoF a homogeneous Poisson process and tho the GoF of a homogeneous Poisson cluster process, and the model with the best fit is returned as the final result. To avoid optimization problems (i.e., obtaining non-realistic parameters for the Poisson cluster models), arguments \emph{lower} and \emph{upper} alllow restricting the range of values that these parameters can attain. If these ranges are not set  by the user, \emph{select.model.gof} will select by default the most extreme and sensible values (e.g., from just one cluster to as many clusters as points in the pattern). As the experience shows that different optimization algorithms provide different results, the basic framework is repeated using both the "L-BFGS-B" and the "Nelder-Mead" algorithms of \link[stats]{optim}. In adition to the initial parameters provided with the argument \emph{parscale}, an aditional otimization using as \emph{parscale} the maximum and minimum values of the \emph{lower} and \emph{upper} parameters range is also tried with the "L-BFGS-B" algorithm.

The goodness of fit is based on the K function (estimated by \link[spatstat]{Kest} or \link[spatstat]{Kinhom}) and is  evaluated following the approach of  Loosmore and Ford (2006) implemented in \link{LF.gof}. This implies computing the sum of squared differences between the observed K function and  the mean of the K functions of \emph{nsim} simulations from the fitted models (the \emph{u} statistic), which means lots of computations and makes the process time comsuming.

Whereas the  model selection approach based on AIC implemented in \link{select.model2} is more appropriated for inference purposes, the approach followed by \emph{select.model.gof} in general is able to select models whose simulations resemble more closely the original pattern.
}
\value{
select.model.gof returns an object of class "\emph{selectedmodgof}",i.e., a list with components:

   \item{gof.u}{vector with the \emph{u} statistic of each model.}
   \item{best.gof }{the minimum \emph{u} value.}
   \item{best.model }{ The best of the fitted models. }
   \item{models }{vector with the names of the fitted models.}
   \item{gof }{A list the results of the \link[ecespa]{LF.gof} test applied to the selected model, i.e., the   \emph{u} value and the associated \emph{p}-value of the GoF test.}
   \item{envelopes }{An \link[spatstat]{envelope} object with envelopes computed form the best model. }
  \item{pp }{ The analyzed point pattern, with the ppp format of spatstat.}
  \item{best.sigma}{The value of the selected bandwith.}

\emph{envelope.selectedmod} returns a an object of class "fv", see \link[spatstat]{fv.object}, which can be printed and plotted directly. Essentially a data frame containing columns 
  \item{r}{the vector of values of the argument r at which the summary function fun has been estimated.}
  \item{obs}{values of the summary function for the data point pattern.}
  \item{lo}{lower envelope of simulations.}
  \item{hi}{upper envelope of simulations.}
  \item{mmean}{estimated theoretical value of the summary function, computed by averaging simulated values.}
 
 \emph{simulate.selectedmod} returns a list of \emph{nsim} point patterns (with the \link[spatstat]{ppp} format of spatstat, simulated according to the best fitted model.
 
  
}
\references{

Chacon-Labella, J., De la Cruz, M. and Escudero, A. Beyond the classical nurse species effect: diversity assembly in a Mediterranean semiarid dwarf  shrubland. \emph{Journal of Vegetaation Science}. Accepted for publication, July 2015.

Loosmore, N.B. and Ford, E.D. (2006) Statistical inference using the G or K point pattern spatial statistics. \emph{Ecology} 87, 1925-1931

}
\author{
Marcelino de la Cruz
}

\examples{
\dontrun{
# Get the data
data(teucrium)

# Define the sequence of r's at which estimate K(r)
  r<-seq(0,1.5, by=0.01)

# Define different standard deviations for the Gaussian kernel 
# to estimate different intensity surfaces
 sigmas <- seq(0.5, 3.5, by=0.25)



# Fit 28 models (1 Poisson, 1 Poisson cluster, 13 inhomogeneous Poisson
#  and 13 inhomogeneous Poisson cluster)  to  teucrium and select the better ones

teucrium.model <- select.model.gof(teucrium, sigmas=sigmas, r=r)

teucrium.model 

# Show the empirical K function, 
#  and the envelopes based on this model

 plot( teucrium.model, sqrt(./pi)-r~r, legend=F, ylab="L(r)", las=1)


# Compute and plot envelopes for the pcf function according to the best fitted model.

teucrium.env <- envelope(teucrium.model, fun=pcf,  nsim=19)
plot(teucrium.env, legend=F)

# simulate 10 point patterns according to the best fitted model

teucrium.simu <- simulate(teucrium.model, nsim=10)

teucrium.simu


}
}

\keyword{ spatial }
