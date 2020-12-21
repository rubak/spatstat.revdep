\encoding{latin1}
\name{K012}
\alias{K012}
\title{ Tests against 'independent labelling' }
\description{
Given a "fixed" point pattern and some process that asign labels (I,J) to another "variable" point pattern,
  \code{K012} estimates the combined bivariate K function between the fixed pattern and every type of the 
variable pattern, and  test that they are independent (i.e. that the labels are randomly assigned,
 irrespectively of the fixed pattern). 
}
\usage{
K012(X, fijo, i, j, nsim = 99, nrank = 1, r = NULL,
	 correction = "isotropic")
}
\arguments{
  \item{X}{ Multitype marked point pattern. An object with the \code{\link[spatstat.geom]{ppp}} format of \pkg{spatstat}.   }
  \item{fijo}{ Number or character string identifying the mark value of the  "fixed" pattern in X }
  \item{i}{ Number or character string identifying the mark value of the  I pattern in X }
  \item{j}{ Number or character string identifying the mark value of the  J pattern in X }
  \item{nsim}{ Number of simulated point patterns to be generated when computing the envelopes.}
  \item{nrank}{ Integer. Rank of the envelope value amongst the \code{nsim} simulated values. 
A rank of 1 means that the minimum and maximum simulated values will be used. }
 \item{r}{ Numeric vector. The values of the argument r at which the K functions  should be evaluated. }
  \item{correction}{ A character item selecting any of the options "border", "bord.modif", "isotropic", "Ripley" or
 "translate". It specifies the edge correction(s) to be applied. }

}
\details{
  This test was developped to answer some questions about the spatial pattern of survival and
mortality of seedlings and its relationships with adult plants in a plant community (De la Cruz \emph{et al. } 2008).
 In order to evaluate the spatial structures of seedlings fates (survive or die), the null hypothesis of random labelling (Cuzick & Edwards 1990, Dixon 2002)
 would be the appropriate one. This kind of pattern is the result of two hierarchical processes: a first one that generates the pattern of points (seedlings)
 and other that assign "labels" (i.e. "die", "survive") to the points. On the other hand,  to analyze the relationships between the spatial pattern of 
emerging seedlings and the pattern of adult plants (two patterns that have been generated independently), independence would be the appropriate
 null hypothesis (Goreaud & Pellisier 2003). However, testing the relationship between the pattern of seedling fates and the pattern of adult plants does
 not completely fit any of the  mentioned hypotheses because, although the pattern of adult plants and the pattern of, e.g., dead seedlings are generated
 independently, their relationship is conditioned by the dependence of the fate "dead" on the locations of emerging seedlings. This implies that 
one can not apply the usual technique of toroidal shifting one pattern over the other to test the independence hypothesis.
 Instead one must permute the label of the focal fate (i.e. survive, die) over the global pattern of seedlings points,
 keeping the locations and labels of adults fixed. This is the method that \code{K012} uses to build the envelopes. The bivariate K functions are
 computed with the Lotwick's and Silverman's (1982) combined estimator (\code{\link{Kmulti.ls}}).
}
\value{
    A list with two elements.
  \item{k01 }{Bivariate K function of the fixed point pattern and the I variable type, with simulation envelopes}
  \item{k02 }{Bivariate K function of the fixed point pattern and the J variable type, with simulation envelopes}

 \item{ }{}

 \item{ }{Each of the above elements is a \code{\link[spatstat.geom]{fv.object}}, essentially a \code{data.frame} with the following items:
}
\item{r }{the values of the argument r at which the functions kave been estimated}
\item{hi }{upper envelope of simulations}
\item{lo }{lower envelope of simulations}
together with the observed corrected estimate of the combined bivariate K function ( \code{iso}, \code{trans}, \code{border}, etc).
}
\references{
Cuzick, J. and Edwards, R. 1990. Spatial clustering for inhomogeneous populations (with discussion). 
 \emph{Journal of the Royal Statistical Society} B \bold{ 52}: 73-104.

 \enc{De la Cruz, M. 2006. Introducción al análisis de datos mapeados o algunas de las (muchas) cosas
 que puedo hacer si tengo coordenadas. \emph{Ecosistemas} 15 (3): 19-39.}{De la Cruz, M. 2006. Introduccion al analisis 
 de datos mapeados o algunas de las (muchas) cosas  que puedo hacer si tengo coordenadas. \emph{Ecosistemas} 15 (3): 19-39.} 

De la Cruz, M., Romao, R.L.,  Escudero, A. and Maestre, F.T. 2008. Where do seedlings go? A spatio-temporal analysis of
 early mortality in a semiarid specialist. \emph{Ecography}, 31(6): 720-730. \url{http://dx.doi.org/10.1111/j.0906-7590.2008.05299.x}.
 
 
Dixon, P. M. 2002. Ripley's K function. In \emph{The encyclopedia of environmetrics} 
(eds. El-Shaarawi, A.H. & Piergorsch, W.W.), pp. 1976-1803. John Wiley & Sons Ltd, NY.

Goreaud, F. and Pelissier, R. 2003. Avoiding misinterpretation of biotic interactions with the intertype K12-function: population independence 
vs. random labelling hypotheses. \emph{J. Veg. Sci.} \bold{14}: 681-692.

Lotwick, H. W. & Silverman, B. W. 1982. Methods for analysing spatial processes of several types of points. 
 \emph{Journal of the Royal Statistical Society} B \bold{44}: 406-413.

}
\author{Marcelino de la Cruz }
\seealso{\code{\link{dixon2002}} for another segregation test, based in the contingency table of counts of nearest neigbors in a marked point pattern.}
\examples{


data(Helianthemum)


# To get 2.5\% and 97.5\% envelopes of 999 random labelings.
## as in fig. 10 of De la Cruz (2006), change nsim= 999 and nrank=25
nsim <-19
nrank <-1



## Test asociation/repulsion between the fixed pattern of adult
## H. squamatum plants and the "variable" pattern of surviving and 
## dead seedlings,


cosa <- K012(Helianthemum, fijo="adultHS", i="deadpl", j="survpl",
             r=seq(0,200,le=201), nsim=nsim, nrank=nrank, correction="isotropic")

plot(cosa$k01, sqrt(./pi)-r~r,  col=c(3, 1, 3), lty=c(3, 1, 3), las=1,
         ylab=expression(L[12]), xlim=c(0, 200), 
	 main="adult HS vs. dead seedlings")

plot(cosa$k02, sqrt(./pi)-r~r, col=c(3, 1, 3), lty=c(3, 1, 3), las=1, 
         ylab=expression(L[12]), xlim=c(0, 200),
 	 main="adult HS vs. surviving seedlings")

}
\keyword{ spatial}
