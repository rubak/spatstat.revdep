\name{locmincon}
\alias{locmincon}
\title{
  Locally Fitted Cluster or Cox Point Process Model
}
\description{
  Fits a Neyman-Scott cluster process or Cox point process model
  using local minimum contrast.
}
\usage{
locmincon(\dots, sigma = NULL, f = 1/4, verbose = TRUE,
                 localstatargs = list(), LocalStats = NULL,
                 tau = NULL)
}
\arguments{
  \item{\dots}{
    Arguments passed to \code{\link[spatstat]{kppm}}
    to determine the template homogeneous model.
  }
  \item{sigma}{
    Standard deviation of Gaussian kernel for local likelihood.
  }
  \item{f}{
    Argument passed to \code{\link[spatstat]{bw.frac}} to
    compute a value for \code{sigma} if it is missing or \code{NULL}.
  }
  \item{verbose}{
    Logical. If \code{TRUE}, print progress reports. 
  }
  \item{localstatargs}{
    Optional. List of arguments to be passed to the local statistic
    \code{\link[spatstat]{localK}},
    \code{\link[spatstat]{localKinhom}},
    \code{\link[spatstat]{localpcf}} or 
    \code{\link[spatstat]{localpcfinhom}}.
  }
  \item{LocalStats}{
    Optional.
    Values of the local statistics, if they have already been computed.
  }
  \item{tau}{
    Optional. Bandwidth for smoothing the fitted cluster parameters.
  }
}
\details{
  The template or homogeneous model is first fitted by
  \code{\link[spatstat]{kppm}}. 
  The statistic used to fit the template model is determined
  (as explained in the help for \code{\link[spatstat]{kppm}})
  by the arguments \code{statistic} and \code{trend}.

  The \emph{local} version of this statistic is then computed.
  If \code{statistic="K"} and \code{trend=~1}
  for example, the template model is fitted
  using the \eqn{K} function \code{\link[spatstat]{Kest}},
  and the local version is the local \eqn{K} function
  \code{\link[spatstat]{localK}}. The possibilities are:

  \tabular{llll}{
    \code{statistic} \tab stationary? \tab
    \emph{template} \tab \emph{local} \cr
    \code{"K"} \tab yes
    \tab \code{\link[spatstat]{Kest}}
    \tab \code{\link[spatstat]{localK}}
    \cr
    \code{"K"} \tab no
    \tab \code{\link[spatstat]{Kinhom}}
    \tab \code{\link[spatstat]{localKinhom}}
    \cr
    \code{"pcf"} \tab yes
    \tab \code{\link[spatstat]{pcf}}
    \tab \code{\link[spatstat]{localpcf}}
    \cr
    \code{"pcf"} \tab no
    \tab \code{\link[spatstat]{pcfinhom}}
    \tab \code{\link[spatstat]{localpcfinhom}}
  }

  These local functions, one for each data point, are then spatially
  averaged, using a Gaussian kernel with standard deviation \code{sigma}.
  Finally the model is fitted to each of the averaged local functions
  to obtain a local fit at each data point.
}
\value{
  Object of class \code{"locmincon"}.
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
  \code{\link{loccit}}
}
\examples{
   X <- redwood[owin(c(0,1), c(-1,-1/2))]
   fit <- locmincon(X, ~1, "Thomas", sigma=0.07)
   fit
}
\keyword{spatial}
\keyword{models}
