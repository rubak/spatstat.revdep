\encoding{latin1}
\name{SF}
\alias{SF}
\alias{SFphylotree}
\alias{SFtraits}
\docType{data}
\title{
San Francisco forest plot.
}
\description{
Point pattern describing the locations of  822 trees (belonging to 113 species) in the San Francisco forest (southern Ecuador). Accompanied by a phylogenetic tree and a file with wood density of some of the species.
}
\usage{
data(SF)
data(SFphylotree)
data(SFtraits)
}
\format{
 \code{SF} is an object of class \code{ppp}  of \pkg{spatstat} representing the point pattern of trees locations, with a data.frame of marks. See \code{\link[spatstat.geom:ppp]{ppp.object}} 
  for details of the format.The dataset has 822 points with the following marks: 
  \describe{
    \item{species}{Species to which each tree belongs}
    \item{dbh}{Diameter at breast height of each tree}
      }
\code{SFphylotree} is a  phylogenetic tree of the class  \link[ape]{phylo} of \pkg{ape}, showing the phylogenetic relationshios among  296 tree species.

\code{SFtraits} is a data.frame with just one column, giving the estimated wood density for some of the 296 tree species.
}
\section{Details}{
This dataset represents the locations of trees with \emph{dbh} >= 5 cm in a polygonal plot of approx. 100 x 70 m in the montane forest near to Reserva \enc{Biológica}{Biologica} San Francisco  (Zamora-Chinchipe, southern Ecuador). These are part of the data collected by Vicuña (2016) and other Ecuadorian botanists and have been  analyzed several times ( e.g., \enc{Chacón}{Chacon} et al.  2014). The coordinates of trees are given in meters and the \emph{dbh} in cm. The phylogenetic tree has been extracted from Phylomatic and calibrated with the BLADJ algorithm of Phylocom. The data about wood density have been compiled from several sources and is expressed in \eqn{g/cm^3}.

}
\source{
\enc{Vicuña, R.  (2016) \emph{Estructura espacial y dinámica del bosque montano del Sur del Ecuador. Interacciones bióticas y limitaciones abióticas}. Tesis Doctoral. Universidad Politécnica de Madrid.}{Vicuna, R.  (2016) \emph{Estructura espacial y dinamica del bosque montano del Sur del Ecuador. Interacciones bioticas y limitaciones abioticas.} Tesis Doctoral. Universidad Politecnica de Madrid.}

}
\references{
\enc{Chacón-Labella et al. (2014) Negative density dependence and environmental heterogeneity effects on tree ferns across succession in a tropical montane forest.
\emph{Perspectives in Plant Ecology, Evolution and Systematics}  16(2): 52-63.}{Chacon-Labella et al. (2014) Negative density dependence and environmental heterogeneity effects on tree ferns across succession in a tropical montane forest. \emph{Perspectives in Plant Ecology, Evolution and Systematics}  16(2): 52-63.}
}

\keyword{datasets}
