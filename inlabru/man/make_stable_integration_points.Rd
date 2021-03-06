% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/integration.R
\name{make_stable_integration_points}
\alias{make_stable_integration_points}
\title{Basic robust integration weights for mesh/polygon intersections}
\usage{
make_stable_integration_points(mesh, bnd, nsub = NULL)
}
\arguments{
\item{mesh}{Mesh on which to integrate}

\item{bnd}{\code{inla.mesh.segment} defining the integration domain}

\item{nsub}{number of subdivision points along each triangle edge, giving
 \code{(nsub + 1)^2} proto-integration points used to compute
the vertex weights
(default \code{9}, giving 100 proto-integration points)}
}
\value{
\code{list} with elements \code{loc} and \code{weight} with
  one integration point for each mesh vertex of triangles overlapping
  the integration domain
}
\description{
Basic robust integration weights for mesh/polygon intersections
}
\author{
Finn Lindgren <\email{finn.lindgren@gmail.com}>
}
\keyword{internal}
