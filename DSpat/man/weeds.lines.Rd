\name{weeds.lines}
\alias{weeds.lines}
\docType{data}
\title{Transect lines from Dubbo weed data}
\description{
Lines sampled in a farming paddock with eight 150m wide transects (75m each side)}
\usage{data(weeds.lines)}
\format{
  A data frame with 8 observations on the following 6 variables.
  \describe{
    \item{\code{label}}{Label of the transect 1 to 8}
    \item{\code{x0}}{x coordinate for the beginning of the line}
    \item{\code{x1}}{x coordinate for the end of the line}
    \item{\code{y0}}{y coordinate for the beginning of the line}
    \item{\code{y1}}{y coordinate for the end of the line}
    \item{\code{width}}{full width of the transect}
  }
}
\details{
These are lines constructed for the \code{\link{weeds}} data that were provided by Melville and Welsh (see reference below) that
were used in the Biometrics paper on distance sampling.

The code used to create the lines was as follows:

\preformatted{
  weeds.lines=data.frame(label=1:8,x0=75+0:7*150,x1=75+0:7*150,
     y0=rep(0,8),y1=rep(1200,8),width=rep(149.9999,8))
}

The line widths were reduced by 0.0001 so the transects do no abut because
\code{spatstat} treats them as overlapping polygons.
}
\references{
Melville, G. J., and A. H. Welsh. 2001. Line transect sampling in small regions. Biometrics 57:1130-1137.
}
\keyword{datasets}
