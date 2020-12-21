\name{Kmulti.ls}
\alias{Kmulti.ls}
\title{ Lotwick's and Silverman's combined estimator of the marked K-function }
\description{
  For a multitype point pattern, calculates the combined estimator of the bivariate \eqn{Kij(r)} and \eqn{Kji(r)} functions.
}
\usage{
Kmulti.ls(X, I, J, r = NULL, corre = "isotropic")
}

\arguments{
  \item{X}{ Multitype marked point pattern. An object with the \code{\link[spatstat.geom]{ppp}} format of \pkg{spatstat}.}
  \item{I}{ Subset index specifying the points of the first pattern.  }
  \item{J}{ Subset index specifying the points of the second pattern. }
  \item{r}{ Numeric vector. The values of the argument r at which the multitype K function \eqn{K^*ij(r)} should be evaluated. }
  \item{corre}{ A character item selecting any of the options "border", "bord.modif", "isotropic", "Ripley" or
 "translate", as described in \code{\link[spatstat.core]{Kest}}. It specifies the edge correction(s) to be applied. }
}
\details{
  As a consequence of edge effects, the estimators  \eqn{Kij(r)} and \eqn{Kji(r)} of the same bivariate pattern could differ.
 \eqn{K^*ij(r)} is the combined estimator defined by Lotwick and Silverman (1982) as \deqn{nj*Kij(r)+ ni*Kji(r) / (ni + nj) ,} 
\eqn{ni} and \eqn{nj} being respectively the number of points in \eqn{I} and \eqn{J}. 
}
\value{
  An object of class "fv"  (see \code{\link[spatstat.geom]{fv.object}}). Essentially a data frame containing numeric columns 
  \item{r }{The values of the argument r at which the function \eqn{K^*ij(r)}  has been estimated}.
  \item{theo }{The theoretical value of \eqn{K*ij(r)} for a marked Poisson process, namely \eqn{pi * r^2}}.

 together with a column or columns named "border", "bord.modif", "iso" and/or "trans", according to the selected edge corrections.
 These columns contain estimates of the function \eqn{K^*ij(r)} obtained by the edge corrections named.

}

\references{ 
Lotwick,H.W. & Silverman, B. W. 1982. Methods for analysing spatial processes of several types of points. 
 \emph{Journal of the Royal Statistical Society} B,  44: 406-413. \url{https://doi.org/10.1111/j.2517-6161.1982.tb01221.x}.
 
 }
\author{ Marcelino de la Cruz }
\note{ \code{Kmulti.ls} is a wrapper for a convenient use of the \code{\link[spatstat.core]{Kmulti}} function of \pkg{spatstat}. 
Please refer to its help page for additional documentation.}
\examples{

data(amacrine)

plot(Kmulti.ls(amacrine, I=amacrine$marks=="on", J=amacrine$marks=="off", 
	 corre="isotropic"), sqrt(./pi)-r~r, main="")

# compare with Kmulti

plot(Kmulti(amacrine, I=amacrine$marks=="on", J=amacrine$marks=="off"),
         sqrt(iso/pi)-r~r, add=TRUE, col=3)

plot(Kmulti(amacrine, J=amacrine$marks=="on", I=amacrine$marks=="off"),
      sqrt(iso/pi)-r~r, add=TRUE, col=4)

}
\keyword{ spatial }
