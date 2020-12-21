\encoding{Latin1}
\name{checktree}
\alias{checktree}
\alias{checktraits}

\title{Check Trait and Phylo Data}
\description{
  Checks consistency of names of species in spatial, phylogenetic and traits data.
}

\usage{
  checktree(tree,  mippp, idar, correct.phylo)
  checktraits(traits, mippp, idar, correct.trait.na, correct.trait)
  
}

\arguments{
 \item{tree}{A phylogenetic tree in \code{phylo} format (\pkg{ape}) or a phylogenetic covariance matrix}
  \item{traits}{A \link[base]{data.frame} of traits, or a distance matrix among species (in \link[stats]{dist} or \link[base]{matrix} format) computed on a data.frame of traits.} 
  \item{mippp}{A multitype (a.k.a. multivariate) marked point pattern . An object with the \link[spatstat]{ppp} format of \pkg{spatstat}, with the names of species as marks.}
  \item{idar}{Character. The name of the idar fucntion to be computed. Either "isar", "ipsvar", "ipsrar", "ipsear", "ipscar", or "imntdar"}
  \item{correct.phylo}{Character. Either \code{"mean"} meaning \emph{"include missing species in the tree with a constant mean phylogenetic covariance"} or \code{"exclude"}, meaning \emph{"exclude missing species in the tree from the analysis"} }
   \item{correct.trait.na}{Logical flag indicating whether \code{NA} values in the matrix of traits should be "corrected": NA values will be assigned the mean trait value.}
   \item{correct.trait}{Character. Either \code{"mean"} or \code{"exclude"}. Species missing in the data.frame of traits will be assigned mean trait values or  will be excluded from the analysis, respectively. }
}

\details{
This functions check for the coincidence of  species in the point pattern and in the phylogenetic or trait data and for the existence of missing data. If \code{correct.phylo="mean"}, species which are absent in the phylogenetic tree or covariance matrix will be included and assigned mean phylogenetic covariance. If \code{correct.phylo="exclude"}, missing species in the tree will be excluded from the analysis (i.e., will not be considered in the computation of the local phylogenetic indices). If \code{correct.trait.na="TRUE"}, NA values for traits in the data.frame of traits will be assigned the mean trait value. If \code{correct.trait="mean"}, missing species in the data.frame of traits will be included and assigned mean trait values.  If \code{correct.trait="exclude"}, missing species in the data.frame of traits will be excluded from the analysis (i.e., will not be considered in the computation of the local functional indices).
}

\section{Warning}{
       The transcription of species names in the multivariate \code{mippp}, in the row names of  the \code{data.frame} of traits (or in the names or \code{dimnames} of  the distance matrices) should be identical. The same applies to the \code{tiplabels} of the phylogenetic tree.
     }


\value{
  \code{checktree} returns a covariance matrix with the appropriate species corrections.   \code{checktraits} returns a data.frame of traits with the appropriate corrections or a Gower distance matrix among the species from the corrected data.frame of traits if   \code{idar="ifdar"}.
}



\author{ Marcelino de la Cruz \email{marcelino.delacruz@urjc.es}}

\keyword{ spatial }
