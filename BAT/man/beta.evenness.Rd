% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BAT.R
\name{beta.evenness}
\alias{beta.evenness}
\title{Beta diversity evenness (Taxon, Phylogenetic or Functional Diversity - TD, PD, FD).}
\usage{
beta.evenness(
  comm,
  tree,
  distance,
  method = "expected",
  func = "camargo",
  abund = TRUE
)
}
\arguments{
\item{comm}{A sites x species matrix, with either abundance or incidence data.}

\item{tree}{An hclust or phylo object representing a phylogenetic or functional tree.}

\item{distance}{A dist or matrix object representing the phylogenetic or functional distance between species. If both tree and distance are missing, taxonomic evenness is calculated.}

\item{method}{Calculate evenness using "expected" values (default) or values based on "contribution" of species to the tree.}

\item{func}{Calculate evenness using "Camargo" (default) or "Bulla" index.}

\item{abund}{A boolean (T/F) indicating whether evenness should be calculated using abundance data.}
}
\value{
Distance matrix between sites.
}
\description{
Difference of evenness between pairs of sites.
}
\details{
This measure is simply the pairwise difference of evenness calculated based on the index of Camargo (1993) or Bulla (1994) using the values of both species abundances and edge lengths in the tree (if PD/FD).

If no tree or distance is provided the result is the original index.
}
\examples{
comm <- matrix(c(1,2,0,0,0,1,1,0,0,0,0,2,2,0,0,1,1,1,1,100), nrow = 4, byrow = TRUE)
distance <- dist(c(1:5), method = "euclidean")
tree <- hclust(distance, method = "average")
beta.evenness(comm)
beta.evenness(comm, tree)
beta.evenness(comm, tree, method = "contribution")
beta.evenness(comm, tree, abund = FALSE)
}
\references{
Bulla, L. (1994) An index of evenness and its associated diversity measure. Oikos, 70: 167-171.

Camargo, J.A. (1993) Must dominance increase with the number of subordinate species in competitive interactions? Journal of Theoretical Biology, 161: 537-542.
}
