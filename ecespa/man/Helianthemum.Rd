\encoding{latin1}
\name{Helianthemum}
\alias{Helianthemum}
\docType{data}
\title{ Spatial point pattern of Helianthemum squamatum adult plants and seedlings }
\description{
  Locations of \emph{H. squamatum} adult plants and seedlings in a 6 m x 7 m plot over gypsum soil 
in \enc{Chinchón}{Chinchon} (near Madrid, Spain). These are part of the data collected by Romao (2003) that have been
 analyzed several times (Escudero \emph{et al.}2005, De la Cruz 2006, De la Cruz et al. \emph{in press.} ).
 The coordinates of the plans are given in cm. 

}
\usage{data(Helianthemum)}
\format{
 An object of class "ppp"  of \code{spatstat} representing the point pattern of plants locations marked by their type. See \code{\link[spatstat]{ppp.object}} 
  for details of the format.The dataset has 866 points with the following levels: 
  \describe{
    \item{adultHS}{adult \emph{H. squamatum} plants}
    \item{deadp}{dying \emph{H. squamatum} seedlings}
    \item{survpl}{surviving \emph{H. squamatum} seedlings}
  }
}
\source{
\enc{Romao, R.L. 2003. \emph{Estructura espacial de comunidades de gipsófitos: interacciones bióticas
 y constricciones abióticas.} Tesis Doctoral. Universidad Politécnica de Madrid.}{Romao, R.L. 2003. \emph{Estructura espacial de comunidades de gipsofitos: interacciones bioticas
 y constricciones abioticas.} Tesis Doctoral. Universidad Politecnica de Madrid.}
}
\references{

  De la Cruz, M. 2006. \enc{Introducción al análisis de datos mapeados o algunas de las (muchas) cosas  que puedo hacer si tengo coordenadas}{Introduccion al analisis  de datos mapeados o algunas de las (muchas) cosas  que puedo hacer si tengo coordenadas}. \emph{Ecosistemas}  15 (3): 19-39. 

De la Cruz, M., Romao, R.L.,  Escudero, A. and Maestre, F.T. 2008. Where do seedlings go? A spatio-temporal analysis of
 early mortality in a semiarid specialist. \emph{Ecography},31(6): 720-730. \url{http://dx.doi.org/10.1111/j.0906-7590.2008.05299.x}.

Escudero, A., Romao, R.L., De la Cruz, M. & Maestre, F. 2005. Spatial pattern and neighbour effects on 
\emph{Helianthemum squamatum} seedlings in a Mediterranean gypsum community. \emph{Journal of  Vegetation  Science},  16: 383-390. \url{https://doi.org/10.1111/j.1654-1103.2005.tb02377.x}.
}
\examples{
data(Helianthemum)
plot(Helianthemum)
}
\keyword{datasets}
