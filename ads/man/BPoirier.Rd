\encoding{latin1}
\name{BPoirier}
\alias{BPoirier}
\docType{data}
\title{Tree spatial pattern in Beau Poirier plot, Haye forest, France}
\description{
Spatial pattern of 162 beeches, 72 oaks and 3 hornbeams in a 1-ha 140 yr-old temperate forest plot in Haye, France. 
}
\usage{data(BPoirier)}
\format{
A list with 8 components:\cr
\code{$rect    } is a vector of coordinates \eqn{(xmin,ymin,xmax,ymax)} of the origin and the opposite corner of a 110 by 90 m rectangular plot.\cr
\code{$tri1    } is a list of vertice coordinates \eqn{(ax,ay,bx,by,cx,cy)} of contiguous triangles covering the denser part of the plot.\cr
\code{$tri2    } is a list of vertice coordinates \eqn{(ax,ay,bx,by,cx,cy)} of contiguous triangles covering the sparser part of the plot.\cr
\code{$poly1   } is a list of vertice coordinates \eqn{(x,y)} of the polygon enclosing \code{BPoirier$tri1}.\cr
\code{$poly2   } is a list of two polygons vertice coordinates \eqn{(x,y)} enclosing \code{BPoirier$tri2}.\cr
\code{$trees   } is a list of tree coordinates \eqn{(x,y)}.\cr
\code{$species } is a factor with 3 levels \eqn{("beech","oak","hornbeam")} corresponding to species names of the trees.\cr
\code{$dbh     } is a vector of tree size (diameter at breast height in cm). 
}
\source{
Pard�, J. 1981. De 1882 � 1976/80 : les places d'exp�rience de sylviculture du hetre en foret domainiale de Haye. \emph{Revue Foresti�re Fran�aise}, 33: 41-64.
}

\references{
Goreaud, F. 2000. \emph{Apports de l'analyse de la structure spatiale en foret temp�r�e � l'�tude et la mod�lisation des peuplements complexes}. Th�se de doctorat, ENGREF, Nancy, France.\cr\cr
P�lissier, R. & Goreaud, F. 2001. A practical approach to the study of spatial structure in simple cases of heterogeneous vegetation. \emph{Journal of Vegetation Science}, 12: 99-108.
}
\examples{
data(BPoirier)
BP.spp <- spp(BPoirier$trees, mark=BPoirier$species, window=BPoirier$rect)
plot(BP.spp)
}
\keyword{datasets}
