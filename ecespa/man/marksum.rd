\name{marksum}
\alias{marksum}
\alias{plot.ecespa.marksum}
\alias{print.ecespa.marksum}
\title{ Mark-sum measure }
\description{
  An exploratory data analysis technique for marked point patterns. The marked point pattern is mapped to a random field for visual inspection.
}
\usage{
marksum(mippp, R = 10, nx = 30, ny = 30)

## S3 method for ploting objects of class 'ecespa.marksum':
\method{plot}{ecespa.marksum}(x, what="normalized",  contour=FALSE, grid=FALSE,
	ribbon=TRUE,col=NULL ,main=NULL,xlab="",ylab="",...)

}
\arguments{
  \item{mippp}{ A marked point pattern. An object with the \code{\link[spatstat]{ppp}} format of \pkg{spatstat}. }
  \item{R}{ Radius. The distance argument \emph{ r} at which the mark-sum measure should be computed }
  \item{nx}{Grid density (for estimation) in the x-side. }
  \item{ny}{ Grid density (for estimation) in the y-side. }
  
  \item{x}{An object of class \code{'ecespa.marksum'}. Usually, the result of applying \code{marksum} to a point pattern.}
  \item{what}{What to plot. One of \code{"marksum"} (raw mark sum measure), \code{"point"} (point sum measure) or \code{"normalized"} (normalized sum measure).} 
   \item{contour}{Logical; if \code{"TRUE"} add contour to map.}
    \item{grid}{Logical; if \code{"TRUE"} add marked grid to map.}
    \item{ribbon}{Logical; if \code{"TRUE"} add legend to map.}
    \item{col}{Color table to use for the map ( see help file on image for details). }
    \item{main}{Text or expression to add as a title to the plot.}
    \item{xlab}{Text or expression to add as a label to axis x.}
    \item{ylab}{Text or expression to add as a label to axis y.}
        
    \item{...}{Additional parameters to \code{\link[spatstat]{Smooth.ppp}}, \code{\link[spatstat]{density.ppp}} or \code{\link[spatstat]{as.mask}}, to control 
                  the  parameters of the smoothing kernel, pixel resolution, etc. }
}
\details{
  Penttinen (2006) defines the \emph{mark-sum measure} as a smoothed summary measuring locally the contribution of points and marks. For any fixed location \eqn{x} within the 
  observational   window and a distance \eqn{R}, the mark-sum measure \eqn{ S[R](x)} equals the sum of the marks of the points within the circle of radius
  \eqn{R} with centre in  \eqn{x}. The \emph{point-sum measure} \eqn{ I[R](x)} is defined by him as the sum of points within the circle of radius \eqn{R} with centre
  in  \eqn{x}, and describes the contribution of points locally near \eqn{x}. The \emph{normalized mark-sum measure} describes the contribution of marks 
  near \eqn{x} and is defined (Penttinen, 2006) as
  \deqn{ S.normalized[R](x) = S[R](x)/I[R](x)}
  This implementation of \code{marksum} estimates the mark-sum and the point-sum measures in a grid of points whose density is defined by \code{nx} and 
  \code{ny}.
}
\value{
  \code{marksum} gives an object of class \code{'ecespa.marksum'}; basically a list with the following elements:
   \item{normalized }{Normalized mark-sum measure estimated in the grid points. } 
  \item{marksum }{Raw mark-sum measure estimated in the grid points. } 
  \item{pointsum }{Point-sum measure estimated in the grid points. } 
  \item{minus}{Point-sum of the grid points. For advanced use only.} 
  \item{grid}{ Grid of points. } 
  \item{nx }{Density of the estimating grid  in the x-side. } 
  \item{ny }{Density of the estimating grid  in the x-side. } 
  \item{dataname }{Name of the ppp object analysed. } 
  \item{R}{ Radius. The distance argument \emph{r} at which the mark-sum measure has been computed. }
  \item{window}{Window of the point pattern.}
  
  \code{plot.ecespa.marksum} plots the selected mark-sum measure.
}
\references{ 
Penttinen, A. 2006. Statistics for Marked Point Patterns. In \emph{The Yearbook of the Finnish Statistical Society}, pp. 70-91. 
}
\author{ Marcelino de la Cruz Rot }

\seealso{ \code{\link{getis}}, related to the point-sum measure, and  \code{\link[spatstat]{markstat}} for designing different implementations. }
\examples{

   
 data(seedlings1)
   
 seed.m <- marksum(seedlings1, R=25)

 # raw mark-sum measure; sigma is bandwith for smoothing
 plot(seed.m, what="marksum", sigma = 5)  

 # point sum measure
 plot(seed.m, what="pointsum", sigma = 5) 
   
 # normalized  mark-sum measure
 plot(seed.m,  what="normalized", dimyx=200, contour=TRUE, sigma = 5) 

# the same with added grid and normalized  mark-sum measure
plot(seed.m,  what="normalized", dimyx=200,
      contour=TRUE, sigma = 5, grid=TRUE)


}  
\keyword{ spatial }
