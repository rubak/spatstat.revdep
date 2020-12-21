\name{unify.owin}
\alias{unify.owin}
\title{Spatial window unit rescaler}
\description{Rescales any \code{\link[spatstat]{owin}} to fall inside the unit square.}
\usage{
unify.owin(W)
}
\arguments{
  \item{W}{
An object of class \code{\link[spatstat]{owin}} giving the spatial window to be transformed.
}
}

\details{
This function is a simple wrapper for \code{\link[spatstat]{affine}} deployed to rescale a supplied \code{\link[spatstat]{owin}} to fall inside the unit square.
}

\value{
The rescaled \code{\link[spatstat]{owin}}.
}

\examples{
W <- Window(chorley)
U <- unify.owin(W)

par(mfrow=c(1,2))
plot(W,axes=TRUE)
plot(U,axes=TRUE)
}
