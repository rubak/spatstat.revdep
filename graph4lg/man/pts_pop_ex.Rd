% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{pts_pop_ex}
\alias{pts_pop_ex}
\title{pts_pop_ex : details on simulated populations}
\format{
An object of class 'data.frame' with the following columns :
\describe{
\item{ID}{Population ID of the 10 populations}
\item{x}{Site longitude (RGF93)}
\item{y}{Site latitude (RGF93)}
}
}
\usage{
pts_pop_ex
}
\description{
Simulation dataset
10 populations located on a simulated landscape
}
\examples{
data("pts_pop_ex")
str(pts_pop_ex)
}
\references{
\insertRef{landguth2010cdpop}{graph4lg}
There are as many rows as there are sampled populations.
}
\keyword{datasets}
