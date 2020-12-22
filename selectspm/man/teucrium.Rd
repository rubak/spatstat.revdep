\name{teucrium}
\alias{teucrium}
\docType{data}
\title{
Teucrium capitatum point pattern
}
\description{
Locations of  \emph{Teucrium capitatum} plants in a dwarf-shrub community in Central Spain.  They are part of a more extensive dataset collected and analysed by Chacon-Labella et al. (2015). The coordinates of the plants are given in meters.
}
\usage{data(teucrium)}
\format{
  An object of class "ppp" representing the point pattern of tree locations. See \link[spatstat.geom:ppp]{ppp.object} for details of the format of a ppp object.
}
\references{
Chacon-Labella, J., De la Cruz, M. and Escudero, A. Beyond the classical nurse species effect: diversity assembly in a Mediterranean semiarid dwarf  shrubland. \emph{Journal of Vegetaation Science}. Accepted for publication, July 2015.
}
\examples{
data(teucrium)
plot(teucrium) 
}
\keyword{datasets}
