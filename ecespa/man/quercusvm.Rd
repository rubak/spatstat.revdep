\encoding{latin1}
\name{quercusvm}
\alias{quercusvm}
\docType{data}
\title{ Alive and dead oak trees }
\description{
  Locations of alive and dead oak trees (\emph{Quercus robur})  in a secondary wood in Urkiola Natural Park (Basque country, north of Spain). 
  This is part of a more extensive dataset collected and analysed by Laskurain (2008). The coordinates of the trees are given in meters. 

}
\usage{data(quercusvm)}
\format{
 An object of class "ppp"  representing the point pattern of tree locations. Entries include
  \describe{
    \item{x}{Cartesian x-coordinate of tree.}
    \item{y}{Cartesian y-coordinate of tree. }
    \item{marks}{factor with two levels indicating the status of each tree (1 = "alive", 0 = "dead").}
  }
 See \code{\link[spatstat]{ppp}}  for details of the format of a ppp object. 
}

\references{
  \enc{Laskurain, N. A. (2008) \emph{Dinámica espacio-temporal de un bosque secundario en el Parque Natural de Urkiola (Bizkaia).}
  Tesis Doctoral. Universidad del País Vasco /Euskal Herriko Unibertsitatea.}{Laskurain, N. A. (2008) \emph{Dinamica espacio-temporal de un bosque secundario en el Parque Natural de Urkiola (Bizkaia).}
  Tesis Doctoral. Universidad del Pais Vasco /Euskal Herriko Unibertsitatea.}
}

\keyword{datasets}
