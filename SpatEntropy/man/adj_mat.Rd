% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adj_matrix.R
\name{adj_mat}
\alias{adj_mat}
\title{Adjacency matrix.}
\usage{
adj_mat(dist.mat, dd0 = 0, dd1)
}
\arguments{
\item{dist.mat}{An upper-triangular matrix of Euclidean distances, as returned by \code{\link[=euclid_dist]{euclid_dist()}}.}

\item{dd0}{Numeric, minimum distance for the neighbourhood/couples/pairs, 0 by default.}

\item{dd1}{Numeric, maximum distance for the neighbourhood/couples/pairs.}
}
\value{
An \eqn{n}x\eqn{n} upper-triangular adjacency matrix
(where \eqn{n} in the data vector length)
with value 1 if two units are neighbours or form a couple/pair, 0 otherwise.
}
\description{
\code{adj_mat} builds an upper-triangular adjacency matrix for the set of points/areas
in a chosen distance range.
}
\details{
The adjacency matrix is a square matrix, with each row corresponding to a point/area, and with
1 values along the row marking the points/areas that are considered 'neighbours' or 'pairing/coupling'.
In the context of spatial entropy, an adjacency matrix may take two roles.
If Karlstrom and Ceccato's entropy is computed, the adjacency matrix identifies
what areas are neighbours, i.e. what areas enter the computation of the
local entropy of a specific area.
If a spatial entropy based on the transformed variable \eqn{Z} is computed (see \code{\link{shannonZ}}
for details on \eqn{Z}), the adjacency matrix
identifies what pairs/couples of points/areas must be considered for the computation, according to
the chosen distance range of interest.
}
\examples{
dist.mat=euclid_dist(cbind(rep(1:5, each=5), rep(1:5,5)))
plot(cbind(rep(1:5, each=5), rep(1:5,5)))
adj_mat(dist.mat, dd1=dist.mat[1,2]) #for the contiguity matrix
adj_mat(dist.mat, 1, 3)

}
