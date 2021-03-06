% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{data_rainforest}
\alias{data_rainforest}
\title{Rainforest tree data.}
\format{A \code{ppp} object (see package \code{spatstat}) with 7251 points, containing:
\describe{
\item{window}{An object of type \code{owin} (see package \code{spatstat}), the 1000x500 metres observation area}
\item{x, y}{Numeric vectors with points' coordinates}
\item{marks}{A character vector matching the tree species to the data points}
}}
\source{
http://www.ctfs.si.edu
}
\usage{
data_rainforest
}
\description{
A marked point pattern dataset about four rainforest tree species.
}
\details{
This dataset documents the presence of tree species over Barro Colorado Island, Panama.
Barro Colorado Island has been the focus of intensive research on lowland tropical
rainforest since 1923 (http://www.ctfs.si.edu). Research identified several tree species
over a rectangular observation window of size 1000x500 metres; the tree species
constitute the point data categorical mark. This dataset presents 4 species with
different spatial configurations: Acalypha diversifolia, Chamguava schippii,
Inga pezizifera and Rinorea sylvatica. The overall dataset has a total number of 7251 points.
The dataset is analyzed with spatial entropy measures in Altieri et al (2017) (references can be
found at \code{SpatEntropy}).
}
\keyword{datasets}
