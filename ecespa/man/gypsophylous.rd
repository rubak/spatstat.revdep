\encoding{latin1}
\name{gypsophylous}
\alias{gypsophylous}
\docType{data}
\title{ Spatial point pattern of a plant community}
\description{
  Locations of  plants  in a gypsophylous plant community in Central Spain. These are part of the data collected by Romao (2003) that have been
 analyzed several times (Escudero \emph{et al.} 2005, De la Cruz 2006).
 The coordinates of the plans are given in cm. 
}
\usage{data(gypsophylous)}
\format{
 An object of class "ppp"  of \code{spatstat} representing the point pattern of plants locations. See \code{\link[spatstat.geom:ppp]{ppp.object}} 
  for details of the format.
}
\source{
\enc{Romao, R.L. 2003. \emph{Estructura espacial de comunidades de gipsófitos: interacciones bióticas
 y constricciones abióticas.} Tesis Doctoral. Universidad Politécnica de Madrid.}{Romao, R.L. 2003. \emph{Estructura espacial de comunidades de gipsofitos: interacciones bioticas
 y constricciones abioticas.} Tesis Doctoral. Universidad Politecnica de Madrid.}
}
\references{
  \enc{De la Cruz, M. 2006. Introducción al análisis de datos mapeados o algunas de las (muchas) cosas
 que puedo hacer si tengo coordenadas}{De la Cruz, M. 2006. Introduccion al analisis 
 de datos mapeados o algunas de las (muchas) cosas  que puedo hacer si tengo coordenadas}. \emph{Ecosistemas}  15 (3): 19-39. 

Escudero, A., Romao, R.L., De la Cruz, M. & Maestre, F. 2005. Spatial pattern and neighbour effects on 
\emph{Helianthemum squamatum} seedlings in a Mediterranean gypsum community. \emph{Journal of  Vegetation  Science},  16: 383-390. \url{https://doi.org/10.1111/j.1654-1103.2005.tb02377.x}.

}
\examples{
\dontrun{

data(gypsophylous)

plot(gypsophylous)

}
}
\keyword{datasets}
