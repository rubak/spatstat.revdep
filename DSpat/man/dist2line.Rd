\name{dist2line}
\alias{dist2line}
\title{ Compute perpendicular distances and projections onto line}
\description{ Calculates perpendicular distances of a point process contained within a strip
 to the center line of the strip they are contained in.  It also computes the positions
 of the objects projected onto the line. This is the inverse of the \code{\link{offset.points}} function.
}
\usage{
dist2line(object.ppp, line.ends)
}
\arguments{
  \item{object.ppp}{ point process for observations in a strip }
  \item{line.ends}{ ends of line x0,y0,x1,y1}
}

\value{
  \item{distVals}{vector of perpendicular distances}
  \item{projection}{dataframe of projected locations on the line}
}
\author{ Devin Johnson }
\seealso{\code{\link{project2line}}}

