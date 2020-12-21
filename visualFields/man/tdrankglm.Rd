\name{tdrankglm}
\alias{tdrankglm}
\title{GLM fit for TD rank curve}
\description{gets a generalized linear model fit for TD rank curve}
\usage{tdrankglm( tdr, familytxt = c( "gaussian" ), link = make.link( "logit" ),
                  rankCentral = NULL, scaleFactor = 52.4 )}
\arguments{
  \item{tdr}{table with TD rank curve}
  \item{familytxt}{family of distributions to use with \code{\link{glm}} fit. See \code{\link{family}}. Default is \code{gaussian}}
  \item{link}{link function to use with \code{\link{glm}} fit. See \code{\link{family}}. Defaulut is \code{logit}}
  \item{rankCentral}{central ranked positions to use in the fit. By default it is \code{NULL}, so that all rank locations are used for the fit}
  \item{scaleFactor}{scale to use for ranks to make them be between 0 and 1. Default is \code{52.4}}
}
\details{details}

\value{
fitted values of the generalized linear model fit for TD rank curve
}

\seealso{ \code{\link{tdrank}}, \code{\link{tdrankperc}} }

\author{Ivan Marin-Franch}