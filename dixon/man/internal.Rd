\name{dixon-internal}
\alias{check}
\alias{ginv}
\alias{mNNinfoc}
\alias{mNNinfo2}
\alias{mNNinfo2b}
\alias{mNNtest}
\alias{NNid}
\title{ Internal functions. }
\description{Internal functions.}
\usage{
check(x, v, l1, l2)
ginv(m)
mNNinfo2(n, R, Q)
mNNinfo2b(n, R, Q)
mNNinfoc(xy, label, nnid = NULL, splancs = TRUE, fortran)
mNNtest(info, obsN = NULL)
NNid(xy, splancs = TRUE)
}
\details{
These fucntions are called by \code{\link{dixon}} and are usually not to be called by the user.
}
\author{ Philip M. Dixon and Marcelino de la Cruz}
\keyword{ internal }
