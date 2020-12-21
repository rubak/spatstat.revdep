\name{methods.locppm}
\alias{methods.locppm}
\alias{as.interact.locppm}
\alias{as.ppm.locppm}
\alias{coef.locppm}
\alias{confint.locppm}
\alias{is.poisson.locppm}
\alias{print.locppm}
\title{Methods for Local Gibbs Models}
\description{
  Methods for various generic functions, for the class
  \code{"locppm"} of locally fitted Gibbs point process models.
}
\usage{
  \method{as.interact}{locppm}(object)

  \method{as.ppm}{locppm}(object)

  \method{coef}{locppm}(object, \dots, which = c("local", "homogeneous"))

  \method{confint}{locppm}(object, parm, level = 0.95, \dots, which = c("local", "homogeneous"))

  \method{is.poisson}{locppm}(x)

  \method{print}{locppm}(x, \dots)
}
\arguments{
  \item{object,x}{
    A locally-fitted  Gibbs point process model (object of class
    \code{"locppm"}). 
  }
  \item{\dots}{
    Additional arguments passed to the default method
    (for \code{confint.locppm}) or ignored (by \code{coef.locppm}).
  }
  \item{which}{
    Character string determining whether to perform calculations
    for the local Gibbs model (\code{which="local"}, the default) or the
    corresponding homogeneous Gibbs model (\code{which="homogeneous"}).
  }
  \item{parm}{
    The parameter or parameters for which a confidence interval is
    desired. A character string or character vector matching the names 
    of \code{coef(object)}, or an index or index vector that can be
    applied to \code{coef(object)}.
  }
  \item{level}{
    Confidence level: a number between 0 and 1.
  }
}
\value{
  \code{as.interact} returns an interaction structure (object of class
  \code{"interact"}). 

  \code{as.ppm} returns a fitted Gibbs model (object of class
  \code{"ppm"}).

  \code{coef} and \code{confint}
  return a numeric vector if \code{which="homogeneous"}
  and an object of class \code{"ssf"} if \code{which="local"}.
  
  \code{is.poisson} returns a logical value.

  \code{print} returns \code{NULL}.
}
\details{
  Objects of class \code{"locppm"} represent locally fitted Gibbs
  point process models.

  The functions documented here provided methods for this class,
  for the generic functions
  \code{\link[spatstat]{as.interact}},
  \code{\link[spatstat]{as.ppm}},
  \code{\link[stats]{coef}},
  \code{\link[stats]{confint}},
  \code{\link[spatstat]{is.poisson}}
  and
  \code{\link{print}}.

  For the \code{coef} and \code{confint} methods, the calculations
  can be performed either on the locally fitted model or
  on its homogeneous equivalent, by changing the argument \code{which}.
}
\references{
  \localpaper

  Baddeley, A., Rubak, E. and Turner, R. (2015)
  \emph{Spatial Point Patterns: Methodology and Applications with R}.
  Chapman and Hall/CRC Press.

}
\seealso{
  \code{\link{locppm}}
}
\examples{
   fit <- locppm(swedishpines, ~1, sigma=9, nd=20,
                 vcalc="full", locations="coarse")
   fit
   is.poisson(fit)
   coef(fit)
   coef(fit, which="homogeneous")
   confint(fit)
   confint(fit, which="homogeneous")
   as.ppm(fit)
   as.interact(fit)
}
\author{Adrian Baddeley}
\keyword{spatial}
\keyword{methods}

