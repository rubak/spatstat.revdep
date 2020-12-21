\encoding{latin1}
\name{seedlings}
\alias{seedlings}
\alias{seedlings1}
\alias{seedlings2}
\docType{data}
\title{ Cohorts of Helianthemum squamatum seedlings }
\description{
  Marked point patterns of two consecutive cohorts of seedlings of \emph{H. squamatum} growing in a gypsophylous plant community in Central Spain. The datasets contains the
  locations of the seedlings marked with their heigth. Both the coordinates and the heigth of the seedlings are given in cm.
}
\usage{
data(seedlings1)
data(seedlings2)
}
\format{
  \code{seedlings1} and \code{seedlings2} are objects of class "ppp" representing the point pattern of seedling locations marked by their heights. See \code{\link[spatstat]{ppp.object}} 
  for details of the format. }

\source{
\enc{Romao, R.L. 2003. \emph{Estructura espacial de comunidades de gipsófitos: interacciones bióticas
 y constricciones abióticas.} Tesis Doctoral. Universidad Politécnica de Madrid.}{Romao, R.L. 2003. \emph{Estructura espacial de comunidades de gipsofitos: interacciones bioticas
 y constricciones abioticas.} Tesis Doctoral. Universidad Politecnica de Madrid.}
}

\references{
 
De la Cruz, M. 2006. \enc{Introducción al análisis de datos mapeados o algunas de las (muchas) cosas  que puedo hacer si tengo coordenadas}{Introduccion al analisis  de datos mapeados o algunas de las (muchas) cosas  que puedo hacer si tengo coordenadas}. \emph{Ecosistemas}  15 (3): 19-39. 

  
Escudero, A., Romao, R.L., De la Cruz, M. & Maestre, F. 2005. Spatial pattern and neighbour effects on 
\emph{Helianthemum squamatum} seedlings in a Mediterranean gypsum community. \emph{Journal of  Vegetation  Science},  16: 383-390. \url{https://doi.org/10.1111/j.1654-1103.2005.tb02377.x}.

}
\examples{
\dontrun{

 data(seedlings1)
 
 plot(seedlings1)

}
}
\keyword{datasets}
