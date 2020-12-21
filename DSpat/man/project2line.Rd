\name{project2line}
\alias{project2line}
\title{ Project points onto line}
\description{
 Projects point process contained in strips to the center line of each strip containing
 points. This is the inverse of the \code{\link{create.points.by.offset}} function.
}
\usage{
project2line(obs.ppp, lines.psp)
}
\arguments{
  \item{obs.ppp}{point process contained in strips}
  \item{lines.psp}{ line segment process with label field }
}
\value{
  dataframe of projected locations (x,y) on the lines
}
\author{ Jeff Laake }
\seealso{\code{\link{dist2line}}}

