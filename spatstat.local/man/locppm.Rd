\name{locppm}
\alias{locppm}
\title{
  Locally Fitted Poisson or Gibbs Point Process Model
}
\description{
  Fits Poisson or Gibbs point process model
  using local likelihood or pseudolikelihood.
}
\usage{
locppm(\dots, sigma = NULL, f = 1/4,
       vcalc = c("none", "t", "hessian", "hom", "lik", "full"),
       locations=c("split", "fine", "coarse"),
       ngrid = NULL, grideps = NULL, verbose = TRUE,
       use.fft=FALSE, fft.algorithm="closepairs")
}
\arguments{
  \item{\dots}{
    Arguments passed to \code{\link[spatstat.core]{ppm}}
    to fit the homogeneous model.
  }
  \item{sigma}{
    Standard deviation of Gaussian kernel for local likelihood.
  }
  \item{f}{
    Argument passed to \code{\link[spatstat.core]{bw.frac}} to
    compute a value for \code{sigma} if it is missing or \code{NULL}.
  }
  \item{vcalc}{
    Type of variance calculation to be performed. See Details.
  }
  \item{locations}{
    Spatial locations for local calculations. See Details.
  }
  \item{ngrid}{
    Dimensions of coarse grid, if used. See Details.
    Incompatible with \code{grideps}.
  }
  \item{grideps}{
    Grid spacing of  coarse grid, if used. See Details.
    Incompatible with \code{ngrid}.
  }
  \item{verbose}{
    Logical. If \code{TRUE}, print progress reports.
  }
  \item{use.fft}{
    Logical value indicating whether to perform
    computations using the Fast Fourier Transform.
    With \code{use.fft = TRUE} the code runs much faster
    but some quantities are not computed exactly.
    See Details.
   }
  \item{fft.algorithm}{
    Developer use only.
  }
}
\details{
  This function fits a Poisson or Gibbs point process model
  to point pattern data by local likelihood or local pseudolikelihood
  respectively.

  This command should be used in the same way as
  \code{\link[spatstat.core]{ppm}}. 
  The point pattern data and the specification of the model
  are given in the leading arguments \code{\dots} which are passed
  directly to \code{\link[spatstat.core]{ppm}}.

  In all cases, the local estimates of the coefficients are
  computed. However, because the variance calculations are
  time-consuming, the default is not to perform them.
  This is controlled by the argument \code{vcalc}.
  
  \describe{
    \item{\code{vcalc = "none"}:}{
      no variance calculations are performed.
    }
    \item{\code{vcalc = "t"}:}{
      the \eqn{t} statistic for each parameter is computed
      for the local model.
    }
    \item{\code{vcalc = "hessian"}:}{
      the local Hessian matrix is computed,
      and its negative inverse is used as a surrogate for the local variance.
    }
    \item{\code{vcalc = "hom"}:}{
      No local fitting is performed.
      Calculations are performed only for the homogeneous (template)
      model. The variance of the local parameter estimates
      \emph{under the homogeneous model} is computed.
    }
    \item{\code{vcalc = "lik"}:}{
      In addition to the calculations for \code{vcalc="hom"} described
      above, if \code{use.fft=FALSE}
      the algorithm also computes the local composite likelihood
      ratio test statistic for the test of homogeneity.
      If \code{use.fft=TRUE} then \code{vcalc="lik"} is equivalent to
      \code{vcalc="hom"}.
    }
    \item{\code{vcalc = "full"}:}{
      all variance calculations are performed for the local model.
    }
  }

  The spatial locations, where the model fits and variance calculations
  are performed, are determined by the argument \code{locations}.

  \describe{
    \item{\code{locations = "fine"}:}{
      The calculations are performed at every quadrature point
      of the model. This can take a very long time.
    }
    \item{\code{locations = "coarse"}:}{
      The calculations are performed at the points of a coarse grid
      with dimensions specified by \code{ngrid} or \code{grideps}.
    }
    \item{\code{locations = "split"}:}{
      The fitted coefficients are computed at every quadrature point
      of the model, but the variance calculations (if any) are
      performed at a coarse grid of locations,
      specified by \code{ngrid} or \code{grideps}.
      If neither \code{ngrid} nor \code{grideps} is specified,
      the default is \code{ngrid=10}.
    }
  }

  If \code{use.fft=FALSE} (the default), all desired quantities are
  computed exactly, by an iterative algorithm that 
  fits a separate model at each spatial location. This can be quite
  slow.
  
  If \code{use.fft=TRUE}, we only compute quantities that can be
  obtained using the Fast Fourier Transform, resulting in much faster
  calculations (sometimes 3 orders of magnitude faster) when
  \code{locations="fine"}.
  Properties of the homogeneous model are
  computed accurately. Properties of the locally-fitted model are
  approximated by a first order Taylor expansion. 
}
\value{
  An object of class \code{"locppm"} representing the fitted model.
}
\references{
  \localpaper

  Baddeley, A., Rubak, E. and Turner, R. (2015)
  \emph{Spatial Point Patterns: Methodology and Applications with R}.
  Chapman and Hall/CRC Press.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{methods.locppm}}, 
  \code{\link{plot.locppm}}
}
\examples{
   fit <- locppm(swedishpines, ~1, sigma=9, nd=20)
   fit
}
\keyword{spatial}
\keyword{models}
