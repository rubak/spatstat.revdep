\name{Allogny}
\encoding{latin1}
\alias{Allogny}
\docType{data}
\title{Spatial pattern of oaks suffering from frost shake in Allogny, France.}
\description{
Spatial pattern of sound and splited oaks (\emph{Quercus petraea}) suffering from frost shake in a 2.35-ha plot in Allogny, France.
}
\usage{data(Allogny)}
\format{
A list with 4 components:\cr
\code{$rect   } is a vector of coordinates \eqn{(xmin,ymin,xmax,ymax)} of the origin and the opposite corner of a 125 by 188 m square plot.\cr
\code{$trees  } is a list of tree coordinates \eqn{(x,y)}.\cr
\code{$status } is a factor with 2 levels \eqn{("splited","sound")}.\cr
}
\source{
  Grandjean, G., Jabiol, B., Bruchiamacchie, M. and Roustan, F. 1990. \emph{Recherche de corr?lations entre les param?tres ?daphiques, et plus sp?cialement texture, hydromorphie et drainage interne, et la r?ponse individuelle des chenes sessiles et p?doncul?s ? la g?livure.} Rapport de recherche ENITEF, Nogent sur Vernisson, France.
}
\references{
Goreaud, F. & P?lissier, R. 2003. Avoiding misinterpretation of biotic interactions with the intertype \emph{K12}-function: population independence vs. random labelling hypotheses. \emph{Journal of Vegetation Science}, 14: 681-692.
}
\examples{
data(Allogny)
allo.spp <- spp(Allogny$trees, mark=Allogny$status, window=Allogny$rect)
plot(allo.spp)
}
\keyword{datasets}
