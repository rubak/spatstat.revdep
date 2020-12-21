\name{haz.ppp}
\alias{haz.ppp}
\title{Easily convert xy data to ppp format}
\description{
  ppp maker for the impatient  layman
}
\usage{
haz.ppp(W)
}
\arguments{
  \item{W}{ a data frame or matrix with two or three columns (coordinate x, coordinate y, and mark of the point)  }
}
\details{
  This naive function easily transform your xy data to the format required by spatstat (version <2.0). 
It establishes the window of observation as the rectangle defined by the xy range.
 It asumes that the first two columns are coordinates x and y, and the third (if any) gives the marks of the points.
}
\value{
A point pattern with the format of spatstat v. < 2.0
}
\note{
This window is an underestimate of the real sampling window of the observed points. For a less biassed estimate, use function \code{\link[spatstat.geom]{ripras}} in spatstat.
}
\author{ Marcelino de la Cruz }
\seealso{ \code{\link[spatstat.geom]{ppp}}, \code{\link[spatstat.geom]{as.ppp}}  }
\examples{


data(fig1)

plot(fig1) #typical xyplot

fig1.ppp <- haz.ppp (fig1)

fig1.ppp

plot(fig1.ppp) # point pattern plot of spatstat


}
\keyword{ spatial }
