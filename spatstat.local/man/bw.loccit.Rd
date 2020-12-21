\name{bw.loccit}
\alias{bw.loccit}
\title{
  Cross Validated Bandwidth Selection for Locally Fitted Point Process Model
}
\description{
  Uses cross-validation to select a smoothing bandwidth
  for locally fitting a Cox or cluster point process model.
}
\usage{
bw.loccit(..., use.fft=TRUE, 
               srange = NULL, ns = 9, sigma = NULL,
               fftopt=list(), 
               verbose = TRUE)
}
\arguments{
  \item{\dots}{
    Arguments passed to \code{\link[spatstat]{kppm}} to fit the homogeneous
    version of the model.
  }
  \item{use.fft}{
    Logical value indicating whether to use a quick-and-dirty
    approximation based on a first order Taylor expansion.
  }
  \item{srange}{
    Range of values of the smoothing parameter \code{sigma}
    to be searched. A numeric vector of length 2 giving the minimum
    and maximum values of \code{sigma}.
  }
  \item{ns}{
    Number of values of the smoothing parameter \code{sigma}
    in the range \code{srange} to be searched. A positive integer.
  }
  \item{sigma}{
    Vector of values of the smoothing parameter to be searched.
  }
  \item{fftopt}{Developer use only.}
  \item{verbose}{
    Logical value indicating whether to display progress reports.
  }
}
\details{
  This function determines the optimal value of the smoothing
  parameter \code{sigma} to be used in a call to \code{\link{loccit}}.

  The function \code{\link{loccit}} fits
  a Cox or cluster point process model
  to point pattern data by local composite likelihood.
  The degree of local smoothing is controlled by a smoothing parameter
  \code{sigma} which is an argument to \code{\link{loccit}}.

  For each value of \code{sigma} in a search interval,
  the function \code{bw.loccit} fits the model locally
  and evaluates a cross-validation criterion. The optimal value of
  \code{sigma} is returned. 
}
\value{
  A numerical value giving the selected bandwidth.
  The result also belongs to the class \code{"bw.optim"}
  which can be plotted.
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
  Ns <- if(interactive()) 16 else 2
  b <- bw.loccit(X, ~1, "Thomas", srange=c(0.07, 0.14), ns=Ns)
  b
  plot(b)
}
\keyword{spatial}
\keyword{models}
\keyword{smooth}
