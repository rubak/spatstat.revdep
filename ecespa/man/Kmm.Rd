\encoding{latin1}
\name{Kmm}
\alias{Kmm}
\alias{plot.ecespa.kmm}
\alias{print.ecespa.kmm}
\alias{ecespa.kmm}
\title{ Mark-weighted K-function }
\description{
  This is a functional data summary for marked point patterns that measures the joint pattern of points and marks at different scales determined by \eqn{r}.
}
\usage{
Kmm(mippp, r = 1:10, nsim=NULL)

## S3 method for ploting objects of class 'ecespa.kmm':
\method{plot}{ecespa.kmm}(x, type="Kmm.n", q=0.025, 
            xlime=NULL, ylime=NULL,  maine=NULL, add=FALSE, kmean=TRUE,
            ylabe=NULL, xlabe=NULL, lty=c(1,2,3), col=c(1,2,3), lwd=c(1,1,1),
             ...)

}
\arguments{
   \item{mippp}{ A marked point pattern. An object with the \code{\link[spatstat.geom]{ppp}} format of \pkg{spatstat}. }
   \item{r}{ Sequence of distances at which \code{Kmm} is estimated. }
   \item{nsim}{ Number of simulated point patterns to be generated when computing the envelopes.}
  
  \item{x}{ An object of class 'ecespa.kmm'. The result of applying \code{Kmm} to a marked point pattern.}
  \item{type}{Type of mark-weighted K-function to plot. One of \code{"Kmm"} ("plain" mark-weighted K-function) or 
                        \code{"Kmm.n"} (normalized mark-weighted K-function). }
  \item{q}{Quantile for selecting the simulation envelopes.}
  \item{xlime}{Max and min coordinates for the x-axis.}
  \item{ylime}{Max and min coordinates for the y-axis.}
  \item{maine}{Title to add to the plot.}
  \item{add}{Logical. Should the kmm.object be added to a previous plot? }
  \item{kmean}{Logical. Should the mean of the simulated Kmm envelopes be ploted?}
  \item{ylabe}{Text  or expression to label the y-axis. }
  \item{xlabe}{Text  or expression to label the x-axis.}
  \item{lty}{Vector with the line type for the estimated Kmm function, the simulated envelopes and the mean of the simulated envelopes. }
  \item{col}{Vector with the color for the estimated Kmm function, the simulated envelopes and the mean of the simulated envelopes.}
  \item{lwd}{Vector with the line width for the estimated Kmm function, the simulated envelopes and the mean of the simulated envelopes.}
  \item{\dots}{Additional graphical parameters passed to plot.}
    
}
\details{
  Penttinnen (2006) defines \eqn{Kmm(r)}, the mark-weighted \eqn{K}-function  of a stationary marked point process \eqn{X}, 
  so that \deqn{lambda*Kmm(r) = Eo[sum(mo*mn)]/mu^2} where \eqn{lambda} is the intensity of the process, 
  i.e. the expected number of points of \eqn{X} per unit area, \eqn{Eo[ ] } denotes expectation (given that there is a point at the origin); 
  \eqn{m0} and \eqn{mn} are the marks attached to every two points of the process separated by a distance \eqn{<= r}  and \eqn{mu}
  is the mean mark. It measures the joint pattern of marks and points at the scales determmined by \eqn{r}. If all the marks are set to 1,
  then \eqn{lambda*Kmm(r)} equals the expected number of additional random points within 
  a distance \eqn{r} of a typical random point of \eqn{X}, i.e. \eqn{Kmm} becomes the conventional Ripley's \eqn{K}-function for unmarked point processes.
  As the \eqn{K}-function measures clustering or regularity among the points regardless of the marks, one can separate clustering of marks with the
  \emph{normalized weighted K-function} \deqn{Kmm.normalized(r) = Kmm(r)/K(r)} If the process is independently marked, \eqn{Kmm(r)} equals \eqn{K(r)} so the 
  normalized mark-weighted \eqn{K}-function will equal 1 for all distances \eqn{r}.
  
  If \code{nsim != NULL}, \code{Kmm}  computes  \emph{'simulation envelopes'} from the simulated point patterns. These are simulated from \code{nsim} random permutations
  of the marks over the points coordinates. This is a kind of pointwise test of \eqn{Kmm(r) == 1 } or \eqn{normalized Kmm(r) == 1} for a given \eqn{r}.  
  }
\value{
\code{Kmm} returns an object of class \code{'ecespa.kmm'}, basically a list with the following items:

\item{dataname }{Name of the analyzed point pattern.}
\item{r }{Sequence of distances at which \code{Kmm} is estimated. }
\item{nsim}{Number of simulations for computing the envelopes, or \code{NULL} if none.}
\item{kmm }{Mark-weighted \eqn{K}-function.}
\item{kmm.n }{Normalized mark-weighted \eqn{K}-function.}
\item{kmmsim }{Matrix of simulated  mark-weighted \eqn{K}-functions, or or \code{NULL} if none.}
\item{kmmsim.n }{Matrix of simulated  normalized mark-weighted \eqn{K}-functions, or or \code{NULL} if none.}  
  
}
\references{
De la Cruz, M. 2008. \enc{Métodos para analizar datos puntuales. 
  En: \emph{Introducción al Análisis Espacial de Datos en Ecología y Ciencias Ambientales: Métodos y Aplicaciones } 
  (eds. Maestre, F. T., Escudero, A. y Bonet, A.), pp 76-127. Asociación Española de Ecología Terrestre, Universidad Rey Juan Carlos
  y Caja de Ahorros del Mediterráneo, Madrid.}{Metodos para analizar datos puntuales. 
  En: \emph{Introduccion al Analisis Espacial de Datos en Ecologia y Ciencias Ambientales: Metodos y Aplicaciones } 
  (eds. Maestre, F. T., Escudero, A. y Bonet, A.), pp 76-127. Asociacion Espanola de Ecologia Terrestre, Universidad Rey Juan Carlos
  y Caja de Ahorros del Mediterraneo, Madrid.}


Penttinen, A. 2006. Statistics for Marked Point Patterns. In \emph{The Yearbook of the Finnish Statistical Society}, pp. 70-91. 
}
\author{ Marcelino de la Cruz Rot }

\note{This implementation estimates \eqn{Kmm(r)} without any correction of border effects, so it must be used with caution. However, as \eqn{K(r)} is 
also estimed without correction it migth compensate the border effects on the normalized \eqn{Kmm}-function.
}
\seealso{\code{\link[spatstat.core]{markcorr}} }
\examples{

  ## Figure 3.10 of De la Cruz (2008):
  # change r to r=1:100
  
   r = seq(1,100, by=5)
  
  data(seedlings1)
  
  data(seedlings2)
  
  s1km <- Kmm(seedlings1, r=r)
  
  s2km <- Kmm(seedlings2, r=r)
  
  plot(s1km, ylime=c(0.6,1.2), lwd=2, maine="", xlabe="r(cm)")

  plot(s2km,  lwd=2, lty=2, add=TRUE )

  abline(h=1, lwd=2, lty=3)
  
  legend(x=60, y=1.2, legend=c("Hs_C1", "Hs_C2", "H0"),
	 lty=c(1, 2, 3), lwd=c(3, 2, 2), bty="n")
\dontrun{
## A pointwise test of normalized Kmm == 1 for seedlings1:

   s1km.test <- Kmm(seedlings1, r=1:100, nsim=99)

   plot(s1km.test,  xlabe="r(cm)")

  }
}
\keyword{spatial }

