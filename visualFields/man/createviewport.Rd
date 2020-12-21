\name{ createviewport }
\alias{ createviewport }
\title{Wrapup for createviewport in package \code{\link{grid}}}
\description{
create regions in where to print graphs. Used for the generation of printouts
}
\usage{createviewport( name, left, top, width, height, pwidth = 8.27,
                pheight = 11.69 )}
\arguments{
  \item{name}{name of he viewport}
  \item{left}{a numeric vector or unit object specifying left location (in inches)}
  \item{top}{a numeric vector or unit object specifying topvflocation (in inches)}
  \item{width}{a numeric vector or unit object specifying width (in inches)}
  \item{height}{a numeric vector or unit object specifying height (in inches)}
  \item{pwidth}{width of the page (in inches). Default is \code{8.27}}
  \item{pheight}{height of the page (in inches). Default is \code{11.69}}
}
\details{
create regions in where to print graphs. Used for the generation of printouts
}

\value{
returns a viewport object
}

\seealso{ \code{\link{vflayout}} }

\author{Ivan Marin-Franch}