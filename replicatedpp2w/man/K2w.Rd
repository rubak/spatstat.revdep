\name{K2w}
\alias{K2w}
\alias{print.k2w}
\alias{plot.k2w}
\title{
Two-Way ANOVA-Like Method to Analyze Replicated Point Patterns
}
\description{
Test for effects of  both individual factors and their interaction on  replicated spatial patterns in a two factorial design.
}
\usage{
K2w(pplist = NULL, dataKijk = NULL, nijk = NULL, r, r0 = NULL, rmax = NULL,
        tratA, tratB = NULL, wt = NULL, nsim = 999, correction = "trans", ...)
\method{print}{k2w}(x,...)
\method{plot}{k2w}(x, trat=NULL, ...,  lty = NULL, col = NULL, 
    lwd = NULL, xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL, 
     legend = TRUE, legendpos = "topleft",  fun="L", main=NULL)
}
\arguments{
  \item{pplist}{
A list of point patterns with the \link[spatstat.geom]{ppp} format of \pkg{spatstat}. This argument os alternative to \code{dataKijk}.
}
  \item{dataKijk}{
A \code{data.frame} with the K-functions (in columns) of each of the replicated point patterns. This argument os alternative to \code{pplist}.
}
  \item{nijk}{
A vector with the number of points in each of the replicated point patterns.
}
  \item{r}{
Vector of values for the argument r at which the K functions have been or should be evaluated. If the K functions are to be computed (i.e., if \code{dataKijk} is NULL), first element of \code{r} should be 0.
}
  \item{r0}{
Minimum  r value for which K(r) estimates will be employed  to compute BTSS.
}
  \item{rmax}{
Maximum  r value for which K(r) estimates will be employed to compute BTSS.
}
  \item{tratA}{
A \code{factor} asigning the levels of the first factor to each point pattern replicate.
}
  \item{tratB}{
A \code{factor} asigning the levels of the second factor to each point pattern replicate.
}
  \item{wt}{
A weighting function employed to compute the BTSS. It can be either an R expression, a vector (with \code{length(wt) ==length (r)}) or a single value (i.e., same weight for all distances). By default, \code{K2w} will use \code{wt=r^-2 }. 
}
  \item{nsim}{
Number of resamples to estimate the sampling distribution of the BTSS statistic.
}
  \item{correction}{
Any selection of the options "border", "bord.modif", "isotropic", "Ripley", "translate", "translation", "none" or "best". It specifies the edge correction to be applied if K functions should be computed.
}
  \item{\dots}{
Additional arguments for \link[spatstat.core]{Kest} function of \pkg{spatstat} (only applies if K functions should be computed) or additional graphical arguments for the \link{matplot} function.
}
\item{x}{
an object of class \code{k2w}.
}
  \item{trat}{
(optional)  Factor employed to compute the averaged K functions that will be ploted. By default, \code{plot.k2w} plots averaged functions for the interaction (i.e., for combination of levels of factors) if two factors have been employed in the analysis or for the levels of the unique factor employed (in the case of one-way analysis). To select  some specific display, use one of these: "tratA" (for the levels of the first factor), "tratB" (for the levels of the second factor)  or  "tratAB" (for the combinations of levels of both factors).
}
  \item{lty}{
  (optional) numeric vector of values of the graphical parameter \code{lty} controlling the line type of each plot. 
}
  \item{col}{
  (optional) numeric vector of values of the graphical parameter \code{col} controlling the colour of each K function.
}
  \item{lwd}{
  (optional) numeric vector of values of the graphical parameter \code{lwd} controlling the line width of each plot. 
}
  \item{xlim}{
  (optional) range of x axis.
}
  \item{ylim}{
  (optional) range of y axis.
}
  \item{xlab}{
  (optional) label for x axis.
}
  \item{ylab}{
  (optional) label for y axis.
}
  \item{legend}{
  Logical flag or \code{NULL}. If \code{legend=TRUE}, the algorithm plots a legend in the top left corner of the plot, explaining the meaning of the different line types and colours. 
}
  \item{legendpos}{
  The position of the legend. Either a character string keyword (see \link{legend} for keyword options) or a pair of coordinates in the format \code{list(x,y)}.
}
  \item{fun}{
  One of  \code{"K"}  or \code{"L"} to select the function to be displayed. By default , L(r) = sqrt(K(r)/pi)-r would be ploted.
}
  \item{main}{
  text to display as the title of the plot. By default, the name of the \code{k2w} object would be shown.
}

}
\details{
This function implements a  extension of  the non-parametric one-way ANOVA-like method of Diggle et al. (1991)  to the two-way case, and particularly to test the effects of the interaction of two factors on the spatial structure of replicated point patterns.  From a set of K functions, it generates weighted averaged K functions for each level and combinations of levels of the factors and computes a statistic analogous to a \emph{between-treatment sum of squares}  (BTSS) in clasical ANOVA.  More details are available in Ramon et al. (in revision).

}
\value{
\code{K2w} returns an object of class \code{k2w}. Basically, a list with components:

  \item{ btss.i }{  Between treatment sum of squares (BTSS) for  factor A.}
  \item{ btss.j }{  BTSS for  factor B.}
  \item{ btss.ij }{  BTSS for the interaction of factors A and B.}
  \item{ btss.i.res }{  Resampled distribution of the BTSS statistic  for  factor A.}
  \item{ btss.j.res }{  Resampled distribution of BTSS for  factor B.}
  \item{ btss.ij.res }{  Resampled distribution of BTSS for the interaction of factors A and B.}
  \item{ KrepA }{ Weighted average of the replicated K functions for each level of factor A.}
  \item{ KrepB }{ Weighted average of the replicated K functions for each level of factor B.}
  \item{ KrepAB }{ Weighted average of the replicated K functions for each combination of levels  of factors A and B.}
  \item{ K0i  }{ Global weighted average (i.e., all K fucntions averaged together).}
  \item{ K0j  }{ Global weighted average (i.e., all K fucntions averaged together).}
  \item{ K0ij  }{ Global weighted average (i.e., all K fucntions averaged together).}
  \item{ Rik }{  Data.frame with the residual functions for factor A.}
  \item{ Rjk }{  Data.frame with the residual functions for factor B.}
  \item{ Rijk }{  Data.frame with the residual functions for the interaction of factors A and B.}
  \item{  nsumA }{ Total number of points among the replicates in each level of factor A.}
  \item{  nsumB }{ Total number of points among the replicates in each level of factor B.}
  \item{  nsumAB }{ Total number of points among the replicates in each combinatipon of levels of factors A and B.}
  \item{ wt }{ Weighting function employed to compute the BTSS.}
  \item{ tratA }{ Factor A.}
  \item{tratB }{ Factor B.}
  \item{ tratAB }{ Factor with the combination of levels of A and B.}
  \item{ dataKijk }{  Data.frame with the empirical, replicated, K-functions.}
  \item{ nijk }{  Vector with the number of points in each replicate.}
  \item{ r }{  Vector of r distances at which  K functions are estimated.}
  \item{ r0 }{ Minimum  r value for which K values have been employed to compute BTSS.}
  \item{ KA.res }{ Data.frame with the weighted average of replicated K functions for each level of factor A, for each of the nsim resamples.}
  \item{ KB.res }{ Data.frame with the weighted average of replicated K functions for each level of factor B, for each of the nsim resamples. }
  \item{ KAB.res }{ Data.frame with the weighted average of replicated K functions for each combination of levels  of factors A and B, for each of the nsim resamples.}
  \item{ nameA }{ name of the R object with factor A.}
  \item{ nameB }{ name of the R object with factor B.}
}

\references{
Diggle, P.J., Nicholas, L. & Benes, F.M. (1991) Analysis of Variance for Replicated Spatial Point Patterns in Clinical Neuroanatomy. \emph{Journal of the American Statistical Association}, 86: 618-625.

\href{http://dx.doi.org/10.1111/ecog.01848}{Ramon, P., De la Cruz, M., Chacon-Labella, J. & Escudero, A. (2016). A new two-way ANOVA-like method for analyzing replicated point patterns in ecology. \emph{Ecography}, 39:1109-1117.}
}
\author{
Marcelino de la Cruz
}

\examples{

# Get the data
data(croton)

croton.2w <- K2w(pplist=croton$list.ppp,  r=seq(0,8, by=0.1),               
               tratA=croton$elevation, tratB=croton$slope, nsim=99)

croton.2w

plot(croton.2w)

plot(croton.2w, "tratB")


}

\keyword{ spatial }
