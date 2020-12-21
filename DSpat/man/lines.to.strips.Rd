\name{lines_to_strips}
\alias{lines_to_strips}
\title{Convert lines to transects (strips) }
\description{Convert lines (center) with transect widths to strips and compute rotation angle
 from vertical. With the intersect function in gpclib, it also now clips the portions
 of the transects that are outside the study area.}
\usage{
lines_to_strips(lines, study.area, width = NULL)
}
\arguments{
  \item{lines}{ dataframe with fields named label,x0,x1,y0,y1 and optionally width }
  \item{study.area}{owin class giving study area window}
  \item{width}{optional; if all lines have the same width it can be specified here}
}
\details{
The function assumes that the intersection of the strip and the study area only
results in a single intersection polygon. That means the entire strip cannot pass
outside the study area and then come back into the study area as in an aerial transect
that passes over water to over land and then back over water.  In this case, the
line should end when it passes out of the sampled area and restarted when back in the
sampled area.
}
\value{
  \item{lines}{a psp class of lines with label and angles for rotation added}
  \item{transects}{a list of dataframes with polygon coordinates}
}
\author{Jeff Laake}

