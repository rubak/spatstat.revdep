\name{with.locmincon}
\alias{with.locmincon}
\alias{with.loccit}
\title{
  Evaluate an Expression for a Locally Fitted Model
}
\description{
  Given a locally-fitted Cox or cluster point process model,
  evaluate an expression involving the fitted cluster parameters.
}
\usage{
  \method{with}{locmincon}(data, ...)

  \method{with}{loccit}(data, ...)
}
\arguments{
  \item{data}{
    An object of class \code{"locmincon"} or \code{"loccit"}
    representing a locally-fitted Cox or cluster point process model.
  }
  \item{\dots}{
    Arguments passed to \code{with.default} specifying the
    expression to be evaluated.
  }
}
\details{
  These are method for the generic function \code{\link{with}} for
  the classes \code{"locmincon"} and \code{"loccit"}.

  An object of class \code{"locmincon"} or \code{"loccit"} represents
  a locally-fitted Cox or cluster point process model.
  It contains a data frame
  which provides estimates of the cluster parameters at each of the
  data points of the original point pattern dataset.

  The expression specified by \code{\dots} will be evaluated
  in this dataframe. If the result of evaluation
  is a data frame with one row for each data point,
  or a numeric vector with one entry for each data point,
  then the result will be an object of class \code{"ssf"}
  containing this information. Otherwise, the result will be
  a numeric vector.
}
\value{
  An object of class \code{"ssf"} or a numeric vector.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{ssf}}
}
\examples{
  example(locmincon)
  with(fit, kappa * sigma2)
  example(locmincon)
  with(fit, kappa * sigma2)
}
\keyword{spatial}
\keyword{models}
\keyword{manip}
\keyword{programming}
