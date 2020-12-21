\name{weeds.covariates}
\alias{weeds.covariates}
\docType{data}
\title{Covariate grid for Dubbo weed data}
\description{
Grid (1 sq meter) of covariates for farm paddock in Dubbo weed data.}
\usage{data(weeds.covariates)}
\format{
  A data frame with 120 x 120 observations on the following 4 variables.
  \describe{
    \item{\code{x}}{x coordinate for mid-point of grid cell}
    \item{\code{y}}{y coordinate for mid-point of grid cell}
    \item{\code{sheep}}{0 if no sheep and 1 if sheep were present on the transect}
    \item{\code{strip}}{transect number 1 to 8}
  }
}
\details{
This is the constructed set of covariates for the farm paddock for the \code{\link{weeds}} data that were provided by Melville and Welsh (see reference below) that
were used in the Biometrics paper on distance sampling.

The code used to create the covariate grid was as follows:

\preformatted{
xx=expand.grid(seq(5,1195,10),seq(5,1195,10))
weeds.covariates=data.frame(x=xx$Var2,y=xx$Var1,
                 sheep=rep(c(0,1),each= 120^2/2),
                 strip=rep(c(1,2,3,4,5,6,7,8),each=120^2/8))
save(weeds.covariates,file="weeds.covariates.rda")
}
}
\references{
Melville, G. J., and A. H. Welsh. 2001. Line transect sampling in small regions. Biometrics 57:1130-1137.
}
\keyword{datasets}
