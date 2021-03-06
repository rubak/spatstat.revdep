% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normmix.R
\name{Count_pts}
\alias{Count_pts}
\title{Counts points in a window}
\usage{
Count_pts(pp, win)
}
\arguments{
\item{pp}{Object of class \code{sppmix} or \code{\link[spatstat.geom]{ppp}}.}

\item{win}{The window of observation as an
object of class \code{\link[spatstat.geom]{owin}}, defining the window of observation.}
}
\value{
An integer representing the number of points from \code{pp} within the window \code{win}.
}
\description{
This function counts the number of points
from a point pattern within a specified window.

For examples see

\url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#Count_pts}
}
\examples{
\donttest{
truemix_surf=rmixsurf(m = 3,lambda=100, xlim=c(-5,5),ylim=c(-5,5))
genPP=rsppmix(truemix_surf)
plotmix_2d(truemix_surf,genPP)
Count_pts(genPP,spatstat.geom::square(1))
Count_pts(genPP,spatstat.geom::square(2))
Count_pts(genPP,spatstat.geom::square(3))
Count_pts(genPP,spatstat.geom::square(4))
Count_pts(genPP,spatstat.geom::square(5))
Count_pts(genPP,spatstat.geom::owin(c(-5,5),c(-5,5)))}

}
\seealso{
\code{\link{rsppmix}},
\code{\link{rmixsurf}},
\code{\link{plotmix_2d}},
\code{\link[spatstat.geom]{square}},
\code{\link[spatstat.geom]{owin}}
}
\author{
Sakis Micheas
}
