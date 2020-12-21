\name{Egig}
\alias{Egig}
% <---------------------------------------------------------------------->
\title{Moments of the Generalized Inverse Gaussian Distribution}
% <---------------------------------------------------------------------->
\description{
  Expected value of X, log(X), 1/X and variance for the generalized
  inverse gaussian distribution. This function has been recycled from the ghyp R package.

}
% <---------------------------------------------------------------------->
\usage{
Egig(lambda, chi, psi, func = c("x", "logx", "1/x", "var"))
}
% <---------------------------------------------------------------------->
\arguments{
  \item{lambda}{A shape and scale and parameter.}
  \item{chi, psi}{Shape and scale parameters. Must be positive.}
  \item{func}{The transformation function when computing the expected value.
             \code{x} is the expected value (default), \code{log x} returns the
             expected value of the logarithm of \code{x}, \code{1/x} returns the
             expected value of the inverse of \code{x} and \code{var} returns the
             variance.}
}
% <---------------------------------------------------------------------->
\details{
   \code{Egig} with \code{func = "log x"} uses
   \code{\link[numDeriv]{grad}} from the \R package \emph{numDeriv}. See
   the package vignette for details regarding the expectation of GIG
   random variables.

}
% <---------------------------------------------------------------------->
\value{
  \code{Egig} gives the expected value
  of either \code{x}, \code{1/x}, \code{log(x)} or the variance if \code{func} equals \code{var}.
}
% <---------------------------------------------------------------------->
\references{

  Dagpunar, J.S. (1989). \emph{An easily implemented generalised inverse
  Gaussian generator.} Commun. Statist. -Simula., \bold{18}, 703--710.

  Michael, J. R, Schucany, W. R, Haas, R, W. (1976).  \emph{Generating
  random variates using transformations with multiple roots}, The
  American Statistican, \bold{30}, 88--90.

}
% <---------------------------------------------------------------------->
\author{David Luethi and Ester Pantaleo}
% <---------------------------------------------------------------------->
\seealso{\code{\link{best.lqr}}}
% <---------------------------------------------------------------------->
\examples{
Egig(lambda = 10, chi = 1, psi = 1, func = "x")
Egig(lambda = 10, chi = 1, psi = 1, func = "var")
Egig(lambda = 10, chi = 1, psi = 1, func = "1/x")
}
% <---------------------------------------------------------------------->
\keyword{distribution}
\keyword{generalized gaussian}