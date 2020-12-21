\name{spatstat.local-package}
\alias{spatstat.local-package}
\alias{spatstat.local}
\docType{package}
\title{
  Local Composite Likelihood
}
\description{
  Extension of the spatstat package, for fitting
  spatial point process models by
  local composite likelihood.
}
\details{
  The main functions are
  \tabular{ll}{
    \code{\link{locppm}} \tab Local likelihood fit of Poisson model \cr
                         \tab Local pseudolikelihood fit of Gibbs model \cr
    \code{\link{locmincon}} \tab Local minimum contrast fit\cr
		          \tab of Neyman-Scott or Cox model \cr
    \code{\link{loccit}} \tab Local composite likelihood fit\cr
                          \tab of Neyman-Scott or Cox model
  }
}
\author{
  \adrian.
}
\references{
  \localpaper

  Baddeley, A., Rubak, E. and Turner, R. (2015)
  \emph{Spatial Point Patterns: Methodology and Applications with R}.
  Chapman and Hall/CRC Press.
}
\keyword{package}
