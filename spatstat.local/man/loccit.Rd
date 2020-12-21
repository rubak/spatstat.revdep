\name{loccit}
\alias{loccit}
\title{
  Locally Fitted Cluster or Cox Point Process Model
}
\description{
  Fits a Neyman-Scott cluster process or Cox point process model
  using a locally-weighted composite likelihood.
}
\usage{
loccit(X, trend = ~1,
       clusters = c("Thomas", "MatClust", "Cauchy", "VarGamma", "LGCP"),
       covariates = NULL,
       \dots,
       diagnostics = FALSE,
       taylor = FALSE,
       sigma = NULL, f = 1/4,
       clustargs = list(), control = list(), 
       rmax,
       covfunargs=NULL, use.gam=FALSE, nd=NULL, eps=NULL,
       niter=3,
       fftopt = list(),
       verbose = TRUE)
}
\arguments{
  \item{X}{
    Point pattern.
  }
  \item{trend}{
    Formula (without a left hand side) specifying the form of the
    logarithm of the intensity.
  }
  \item{clusters}{
    Character string determining the cluster model.
    Partially matched.
  }
  \item{covariates}{
    The values of any spatial covariates (other than the Cartesian
    coordinates) required by the model.
    A named list of pixel images, functions, windows or numeric constants.
  }
  \item{diagnostics}{
    Whether to perform auxiliary calculations
    in addition to the local estimates of the model parameters.
  }
  \item{\dots}{
    Additional arguments passed to
    \code{\link[spatstat]{as.mask}} to control the spatial resolution
    in the Fast Fourier Transform.
  }
  \item{taylor}{
    Logical value indicating whether to fit the model
    exactly at each spatial location (\code{taylor=FALSE}, the default)
    or to compute a first-order Taylor approximation to the
    fitted parameters (\code{taylor=TRUE}). The Taylor approximation
    is much faster.
  }
  \item{sigma}{
    Standard deviation of Gaussian kernel for local likelihood.
  }
  \item{f}{
    Argument passed to \code{\link[spatstat]{bw.frac}} to
    compute a value for \code{sigma} if it is missing or \code{NULL}.
  }
  \item{clustargs}{
    List of additional parameters for the cluster model,
    passed to the function \code{RFcov} in the
    \code{RandomFields} package.
  }
  \item{control}{
    List of control arguments passed to the generic optimisation
    function \code{\link[stats]{optim}}.
  }
  \item{rmax}{
    Maximum distance between pairs of points that will contribute
    to the composite likelihood.
  }
  \item{covfunargs,use.gam,nd,eps}{
    Arguments passed to \code{\link[spatstat]{ppm}}
    to control the intensity model and intensity fitting.
  }
  \item{niter}{
    Number of iterations in algorithm if \code{taylor=FALSE}.
  }
  \item{fftopt}{Developer use only.}
  \item{verbose}{
    Logical. If \code{TRUE}, print progress reports. 
  }
}
\details{
  This function fits a Cox or cluster process model to point pattern
  data locally, using the local Palm likelihood technique
  (Baddeley, 2016, section 8).

  It can be used in the same way as \code{\link[spatstat]{kppm}}
  and effectively performs local fitting of the same model.
}
\value{
  An object of class \code{"loccit"}.
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
  \code{\link{locppm}}
}
\examples{
   X <- redwood[owin(c(0,1), c(-1,-1/2))]
   fit <- loccit(X, ~1, "Thomas", nd=5, control=list(maxit=20))
   fit
}
\keyword{spatial}
\keyword{models}
