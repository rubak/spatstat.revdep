\name{create.lines}
\alias{create.lines}
\title{Create a systematic sample of parallel lines across a grid}
\description{
Create a systematic set of lines to sample a rectangular grid. The grid is positioned
with a random start on the study area.  The systematic
grid can be set at any angle and the number of lines is set by the spacing
or the spacing is set by width and number of lines.  This is a wrapper function
for \code{rlinegrid} in \code{spatstat}.
}
\usage{
create.lines(study.area,nlines=NULL,width,spacing=NULL,angle=0)
}
\arguments{
  \item{study.area}{owin class defining area}
  \item{nlines}{number of lines}
  \item{width}{full transect width}
  \item{spacing}{spacing distance between centerlines}
  \item{angle}{angle of rotation in degrees anticlockwise from x-axis}
}
\value{lines dataframe with label,x0,y0,x1,y1,width where x0,y0 is beginning
       and x1,y1 is end of the line
}
\author{Jeff Laake}
\seealso{ \code{\link{simCovariates}},\code{\link{simPts}}}
\examples{
study.area=owin(xrange=c(0,100),yrange=c(0,100))
xp=create.lines(study.area,nlines=10,width=5,angle=180)
ls=lines_to_strips(xp,study.area)
plot(ls$lines,lty=2)
plot(owin(poly=ls$transects),add=TRUE)
}

