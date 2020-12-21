\encoding{latin1}
\name{syrjala}
\alias{syrjala}
\alias{syrjala0}
\alias{syrjala.test}
\alias{plot.syrjala.test}
\alias{print.syrjala.test}
\alias{plot.ecespa.syrjala}
\alias{print.ecespa.syrjala}


\title{ Syrjala's test for the difference between the spatial distributions of two populations }
\description{
  Computes a two-sample Cramer-von Mises (and Kolmogorov-Smirnov) type test for a difference between the spatial distributions of two populations.
  It is designed to be sensitive to differences in the way the populations are distributed across the study area but insensitive to differences in abundance between the two populations.
}
\usage{
syrjala0(coords, var1, var2, nsim, R=FALSE)
syrjala(coords = NULL, var1 = NULL, var2 = NULL, nperm = 999)
syrjala.test(ppp1, ppp2, nsim = 999)
\method{plot}{syrjala.test}(x, coline=1, ...)
\method{plot}{ecespa.syrjala}(x, ...)
}

\arguments{
  \item{coords}{A \code{data.frame} with `$x` and `$y` components.}
  \item{var1}{The first numeric variable }
  \item{var2}{The second numeric variable.}
  \item{nperm}{Number of permutations.}
  \item{nsim}{Number of permutations.}
  \item{R}{Logical. Should be computed using R approach?}
   \item{ppp1}{ A marked point pattern, with the \code{\link[spatstat]{ppp}} format of spatstat, 
                      representing the values of some parameter measured on the corresponding sampling locations. }
  \item{ppp2}{ A marked point pattern, with the \code{\link[spatstat]{ppp}} format of spatstat, 
                     representing the values of some other parameter measured on the same locations than \code{ppp1}. }
  \item{x}{An object of class '\code{syrjala.test}' or  '\code{ecespa.syrjala}' resulting from \code{syrjala} or \code{syrjala.test}, respectively.}
  \item{coline}{color for drawing the statistic's line in the plot.}
 \item{...}{ Graphical parameters passed to \code{\link{hist}}.}
}
\details{
  The null hypothesis of  Syrjala's test is that across the study area, the normalized distributions of the two populations are the same (Syrjala, 1996). Population density data are collected at \eqn{K} sampling locations on two populations. 
  Let \eqn{(xk, yk)} denote the coordinates of the \eqn{kth} sampling location (\eqn{ k= 1,...,K} ); 
  let \eqn{d.i(xk, yk)} denote the sample density at the \eqn{Kth} sampling location of the \eqn{ith} population.
   To construct a test that is independent of the population sizes, the observed density data is first normalized:
   \deqn{gamma.i(xk, yk) = di(xk, yk)/Di,} where \eqn{Di} is the sum  of \eqn{d.i(xk, yk)} observations  across the \eqn{K} sampling locations. 
   The value of the cumulative distribution function at the location \eqn{(xk, yk)}  for the \eqn{ith} population, denoted \eqn{GAMMA.i(xk, yk)},
   is the sum of all normalized density observations, \eqn{gamma.i(xk, yk)}, whose location \eqn{(x, y)} is such that 
   \eqn{x <= xk} and \eqn{y <= yk}. The statistic proposed by Syrjala to test the null hypothesis is the square of the difference 
   between the cumulative distribution functions  of the two populations,  summed over all sampling locations,
   that is \deqn{psi = sum{GAMMA.1(xk, yk) - GAMMA.2(xk, yk)}^2.} 
   As \eqn{psi} is not invariant with respect to the 'corner'  of the rectangle enclosing the study area that is chosen as the origin of the coordinate sytem,
   \eqn{psi}  is computed  four times, one with each corner as the origin, and the average \eqn{psi} is employed as the test statistic.
   The level of significance of the observed \eqn{psi} is determined from its position  in the ordered set of test statistic values from all \eqn{2^K} pairwise 
   permutations (that is approximated from a large number of randomly selected permutations).
}
\value{
  Functions \code{syrjala} or \code{syrjala0} (with the argument \code{R=FALSE}) return an object of class '\code{syrjala.test}'. 
  Functions \code{syrjala.test} or \code{syrjala0} (with the argument \code{R=TRUE}) return an object of class '\code{ecespa.syrjala}'.
  In Both cases, the result is a  list with the following elements: 
  \item{cvm.obs}{(class syrjala.test). The observed (averaged) \eqn{psi} statistic for the CvM test.}
  \item{cvm.sim}{(class syrjala.test). A numeric vector with the \code{nperm+1} simulated \eqn{psi}'s statistics (including \code{cvm.obs}).}
  \item{ks.obs}{(class syrjala.test). The observed (averaged) \eqn{psi} statistic for the K-S test.}
  \item{ks.sim}{(class syrjala.test). A numeric vector with the \code{nperm+1} simulated \eqn{psi}'s statistics (including \code{ks.obs}).}
  \item{datanames}{(class syrjala.test). A character vector with the names of the two patterns, the spatial congruence of which is been analyzed.}
  \item{nperm }{(class syrjala.test). The number of permutations employed in the test (not counting the original data).}
  \item{psi.obs }{(class ecespa.syrjala).The observed (averaged) \eqn{psi} statistic.}
  \item{psi.sim }{(class ecespa.syrjala). A vector with the \code{nsim} simulated \eqn{psi}'s statistics.}
  \item{datanames }{(class ecespa.syrjala). A vector with the names of the two point patterns whose spatial congruence is been analyzed.}
  \item{nsim }{(class ecespa.syrjala). The number of permutations employed in the test.}
  
  Both S3 plot methods plot an histogram with the distribution of the simulated \eqn{psi}'s statistics and draws the observed \eqn{psi} as a vertical line.
  
}


\note{  \code{syrjala} or \code{syrjala0} (with the argument \code{R=FALSE}) implement a Fortran version of Syrjala's test. They run considerably faster
than the "whole-R" implementation of \code{syrjala.test} or \code{syrjala0} (with the argument \code{R=TRUE}). This last implementation is supplied for illustrative purposes and to 
maintain compability with previous versions of package \code{ecespa}. One can use function \code{\link{haz.ppp}} to easily build the \code{ppp} objects from a \code{data.frame} with only three columns (x-coordinate,
y-coordinate and abundance).  

This function has been employed to compute Syrjala's test in Rey-Benayas et al. (2008).
}
\section{Warning }{The test requires both populations being sampled in exactly the same sampling locations. 
Althoug this implementation employs \code{\link[spatstat]{ppp}}'s as the supporting data format, this kind of data are \bold{not} spatial point patterns.
They cannot be analysed with the usual tools employed for marked point patterns.} 

\references{ 

 \enc{Rey-Benayas, J.M., de la Montaña, E., Pérez-Camacho, L., de la Cruz, M., Moreno, D., Parejo, J.L. and Suárez-Seoane, S.}{Rey-Benayas, J.M., de la Montana, E., Perez-Camacho, L., de la Cruz, M., Moreno, D., Parejo, J.L. and Suarez-Seoane, S.} 2010. 
Inter-annual dynamics and spatial congruence of a nocturnal bird assemblage inhabiting a Mediterranean agricultural mosaic. \emph{Ardeola},57(2): 303-320.

Syrjala, S. E. 1996. A statistical test for a difference between the spatial distributions of two populations. \emph{Ecology} 77: 75-80. \url{ https://doi.org/10.2307/2265656}.

}
\author{ Jose M. Blanco-Moreno \email{jmblanco@ub.edu} for the Fortran implementation of Syrjala's original QBasic function, Marcelino de la Cruz Rot for the R version, the wrapping functions and the documentation}

\examples{

   
   data(syr1); data(syr2); data(syr3)
   
   coords <- data.frame(x=syr1$x, y=syr1$y); var1<- syr1$marks; var2<- syr2$marks
   
   (stest<-syrjala(coords, var1, var2, 9999))
   plot(stest)
   
   
   \dontrun{
   
   plot(syrjala.test(syr1, syr2, nsim=999)) 
   
   plot(syrjala.test(syr1, syr3, nsim=999)) 
   
   syrjala0(coords, var1, var2, 9999)
   
   syrjala0(coords, var1, var2, 999, R=TRUE)
   
   
   coords <- expand.grid(x=1:10,y=1:10)
   var1 <- runif(100)
   var2 <- runif(100)
   syrjala(coords, var1, var2, 9999)
   
   
}
}

\keyword{ spatial }
\keyword{ distribution }
