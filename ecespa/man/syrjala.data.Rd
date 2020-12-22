\name{syrjala.data}
\alias{syr1}
\alias{syr2}
\alias{syr3}
\docType{data}
\title{ Syrjala test data}
\description{
  Artificial data to exemplify Syrjala's test.
}
\usage{
data(syr1)
data(syr2)
data(syr3)}
\format{
 \code{syr1}, \code{syr2} and \code{syr3} are marked point patterns 
 of class "ppp" representing the coordinates of some sampling locations, 'marked' by 
 the value of some parameters (e.g. density of individuals) measured on them. See \code{\link[spatstat.geom:ppp]{ppp.object}} 
  for details of the format. On the other hand, one can use function \code{\link{haz.ppp}} to easily build \code{ppp} objects 
  appropriate for use with \code{\link{syrjala.test}}.
}

\examples{
data(syr1)
}
\keyword{datasets}
