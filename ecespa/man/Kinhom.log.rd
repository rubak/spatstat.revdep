\name{Kinhom.log}
\alias{Kinhom.log}

\title{ Simulation envelopes from the fitted values of a logistic model }
\description{
  Computes simulation envelopes for (in-)homogeneous K-function simulating from a vector of probabilitiesn.
}
\usage{
Kinhom.log (A, lambda=NULL, mod=NULL, lifemark="0", prob=NULL,
			r=NULL, nsim=99, correction="trans", ngrid=200)


}

\arguments{
  \item{A}{ A marked point pattern with the \code{\link[spatstat.geom]{ppp}} format of \code{spatstat}. }
\item{lambda}{ Optional. Values of the estimated intensity function as a pixel image (object of class
"\code{\link[spatstat.geom]{im}}" of spatstat) giving the intensity values at all locations of \code{A}. } 
 \item{mod}{ A fitted model. An object of class \code{\link[spatstat.core]{ppm}}. }
  \item{lifemark}{ Level of  the marks of  \code{A} which represents the "live" or "succes" cases.}
  \item{prob}{Numeric vector, with length equal to the number of points of \code{A}, represeting 
  the fitted values of a logistic model fitted to \code{A} marks.}
  \item{r}{ Numeric vector. The values of the argument \eqn{r} at which the \eqn{K(r)} functions  should be evaluated. }
  \item{nsim}{ Number of simulated point patterns to be generated when computing the envelopes. }
   \item{correction}{ A character item selecting any of the options "border", "bord.modif",  or "translate". It specifies 
  the edge correction to be applied when computing K-functions. }
    \item{ngrid}{ Dimensions (ngrid by ngrid) of a rectangular grid of locations where \code{\link[spatstat.core]{predict.ppm}} 
  would evaluate the spatial trend of the fitted models. }
   
}
\details{
  
This function is a wrapper to compute the critical envelopes for Monte Carlo test of goodness-of-fit of (in-)homogeneous K functions, 
simulating from the fittted values of a logistic model (i.e. a binomial GLM with logit link) fitted to the marks ("failure", "success") of a "binomially"-marked 
point pattern. This is particularly interesting in plant ecology when considering alternatives to the \emph{random mortality hypothesis} (Kenkel 1988).
This hypothesis is usually tested building Monte Carlo envelopes from the "succesful" patterns resulting from a  random labelling of a "binomially"-marked point pattern 
(this is equivalent to a random thinning of the whole pattern irrespective of  the marks). 
As tree mortality is rarely random but instead can be modelled as a function of a certain number of covariates, the most natural alternative to 
the \emph{random mortality hypothesis} is the \emph{logistic} mortality hypothesis, that can be tested thinning the original pattern of trees
with retention probabilities defined by the fitted values of a logistic model (Batista and Maguire 1998, Olano et al. 2008).

\code{Kinhom.log} will compute the envelopes by thinning the unmarked point pattern \code{A} with retention probabilities \code{prob}. If no 
\code{prob} vector is provided, all points will be thinned with the same probability ( number of "live" points / number of points ), i.e. \code{Kinhom.log}
will compute random thinning envelopes.

\code{Kinhom.log} will compute envelopes both to homogeneous and inhomogeneous K functions. If no \code{lambda} or \code{mode} arguments
are provided, \code{Kinhom.log} assumes that the original pattern is homogeneous and will use a constant \code{lambda} to compute the inhomogeneous K 
(i.e. it will compute the homogeneous K). The most convenient use with inhomogeneous point patterns is to provide the argument \code{mod}
with an inhomogeneous Poisson model fitted to the original pattern of 'live' points (with spatstat function \code{\link[spatstat.core]{ppm}}; see the examples).
This model will be used to compute (and to update in the simulations) the inhomogeneous trend (i.e. the "lambda") of the patterns.  
If the argument \code{lambda} is provided but not \code{mod},  these lambda will be used as a covariate to fit an inhomogeneous Poisson model that 
will be used to compute (and to update in the simulations) the inhomogeneous spatial trend.

\code{Kinhom.log} will produce an object of class 'ecespa.kci' that can be easily ploted (see the examples). 
This is accomplished by  the S3 ploth method \code{\link{plot.ecespa.kci}}; it will plot the K-function and its envelopes
(actually, it will plot the most usual L-function = \eqn{sqrt[K(r)/pi]-r}).

}
\value{
  \code{Kinhom.log} returns an object of class \code{\link{ecespa.kci}}, basically a list with the following items:
   \item{r }{Numeric vector. The values of the argument \eqn{r} at which the \eqn{K(r)} functions  have been evaluated.}
  \item{kia }{Numeric vector. Observed (in-)homogeneous K function.}
   \item{kia.s }{ Matrix of simulated (in-)homogeneous K functions.}
  \item{datanamea }{ Name of  point pattern \code{A}.}
  \item{modnamea }{ Name of model \code{mod}.}
  \item{type }{ Type of analysis. Always "Kinhom.log".}
  \item{probname }{ Name of the vector of fitted retention probabilities \code{prob}.}
  \item{modtrend }{ Spatial trend (formula) of the model \code{mod}.}
  \item{nsim}{ Number of simulations. }
  
  
}
\references{

Batista, J.L.F. and Maguire, D.A. 1998. Modelling the spatial structure of tropical forests. \emph{For. Ecol. Manag.},  110: 293-314. 

Kenkel, N.C. 1988. Pattern of self-thinning in Jack Pine: testing the random mortality hypothesis. \emph{Ecology},  69: 1017-1024.
 
Olano, J.M., Laskurain, N.A., Escudero, A. and De la Cruz, M. 2009. Why and where adult trees die in a secondary temperate forest? 
The role of neighbourhood. \emph{Annals of Forest Science}, 66: 105. \url{http://dx.doi.org/10.1051/forest:2008074}.


}

\author{ Marcelino de la Cruz Rot }

\section{Warning }{As this implementation involves the use of images as the means of evaluation of the (inhomogeneous) spatial trend, and a mask based on
those images will be used as the point pattern window, the "Ripley's" or "isotropic" edge correction can not be employed.

}


\examples{
  
   
   data(quercusvm)
   
   # set the number of simulations (nsim=199 or larger for real analyses)
   nsim<- 19

   # read fitted values from logistic model:
   
   
   probquercus <-c(0.99955463, 0.96563477, 0.97577094, 0.97327199, 0.92437309,
   0.84023396, 0.94926682, 0.89687281, 0.99377915, 0.74157478, 0.95491518,
   0.72366493, 0.66771787, 0.77330148, 0.67569082, 0.9874892, 0.7918891, 
   0.73246803, 0.81614635, 0.66446411, 0.80077908, 0.98290508, 0.54641754,
   0.53546689, 0.73273626, 0.7347013, 0.65559655, 0.89481468, 0.63946334,
   0.62101995, 0.78996371, 0.93179582, 0.80160346, 0.82204428, 0.90050059,
   0.83810669, 0.92153079, 0.47872421, 0.24697004, 0.50680935, 0.6297911, 
   0.46374812, 0.65672284, 0.87951682, 0.35818237, 0.50932432, 0.92293014,
   0.48580241, 0.49692053, 0.52290553, 0.7317549, 0.32445982, 0.30300865,
   0.73599359, 0.6206056, 0.85777043, 0.65758613, 0.50100406, 0.31340849, 
   0.22289286, 0.40002879, 0.29567678, 0.56917817, 0.56866864, 0.27718552,
   0.4910667, 0.47394411, 0.40543788, 0.29571349, 0.30436276, 0.47859015,
   0.31754526, 0.42131675, 0.37468782, 0.73271225, 0.26786274, 0.59506388, 
   0.54801851, 0.38983575, 0.64896835, 0.37282031, 0.67624306, 0.29429766,
   0.29197755, 0.2247629, 0.40697843, 0.17022391, 0.26528042, 0.24373722,
   0.26936163, 0.13052254, 0.19958585, 0.18659692, 0.36686678, 0.47263005,
   0.39557661, 0.68048997, 0.74878567, 0.88352322, 0.93851375)
   
  

   ################################ 
   ## Envelopes for an homogeneous point pattern:
   
   cosap <- Kinhom.log(A=quercusvm, lifemark="0",  prob=probquercus, nsim=nsim)

   plot(cosap)

   
   ################################ 
   ## Envelopes for an inhomogeneous point pattern:
   
   ## First, fit an inhomogeneous Poisson model to alive trees :
   
   quercusalive <- unmark(quercusvm[quercusvm$marks == 0])

    mod2 <- ppm(quercusalive, ~polynom(x,y,2))

    ## Now use mod2 to estimate lambda for K.inhom:
    
    cosapm <- Kinhom.log(A=quercusvm, lifemark="0", prob=probquercus, 
                                   nsim=nsim, mod=mod2)

   
   ################################ 
   ## An example of homogeneous random thinning:
      
   cosa <- Kinhom.log(A=quercusvm, lifemark="0", nsim=nsim)
   
   plot(cosa)
    
    

}

\keyword{ spatial }


