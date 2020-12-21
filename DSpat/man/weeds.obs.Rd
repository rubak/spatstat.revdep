\name{weeds.obs}
\alias{weeds.obs}
\docType{data}
\title{Observations from Dubbo weed data}
\description{
Observed devils claw in a farming paddock from eight 150m wide transects (75m each side).
These are the records from \code{\link{weeds}} that were seen.}
\usage{data(weeds.obs)}
\format{
  A data frame with 479 observations on the following 4 variables.
  \describe{
    \item{\code{label}}{Label of the transect 1 to 8}
    \item{\code{x}}{x coordinate in the farming paddock}
    \item{\code{y}}{y coordinate created randomly for the data}
    \item{\code{distance}}{perpendicular distance from line to weed}
  }
}
\details{
These are the data constructed from \code{\link{weeds}} that were provided by Melville and Welsh (see reference below) that
were used in the Biometrics paper on distance sampling.

The code used to create the data from weeds was as follows:

\preformatted{
data(weeds.all)
weeds.obs=weeds.all[weeds.all$Seen==1,]
weeds.obs$Seen=NULL
save(weeds.obs,file="weeds.obs.rda")
}

}
\references{
Melville, G. J., and A. H. Welsh. 2001. Line transect sampling in small regions. Biometrics 57:1130-1137.
}
\keyword{datasets}
