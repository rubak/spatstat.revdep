\name{offset.points}
\alias{offset.points}
\title{Offset points from the line to actual position}
\description{Convert x,y point locations on the line and a distance (negative is left of
 line for the direction of travel) for a series of points in a strip.}
\usage{
offset.points(line, pts)
}
\arguments{
  \item{line}{vector with named components of x0,y0,x1,y1; line traverses from (x0,y0) to (x1,y1)}
  \item{pts}{dataframe of x,y,distance for each observed  point; x,y is the location
        on the line that is perpendicular to the object; a negative distance
        implies it is on the left side of the line as defined by the direction of
        travel}
}
\value{
pts dataframe with x,y locations of the objects offset from the line at
        the appropriate distance and side.}
\author{Jeff Laake}
\seealso{ \code{\link{create.points.by.offset}}}

