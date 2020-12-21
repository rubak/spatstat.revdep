\name{stimLoc}
\alias{stimLoc}
\title{stimulus location and size}
\description{plots the stimulus location and size}
\usage{stimLoc( perimetry, pattern, eye, txtfont = "mono", pointsize = 7,
         xminmax = 29, yminmax = 29 )}
\arguments{
  \item{perimetry}{type of perimetry}
  \item{pattern}{pattern used in the perimetry}
  \item{eye}{eye being tested}  
  \item{txtfont}{font of the text with visual-sensitivity values. Default is \code{mono}}
  \item{pointsize}{size of the text with visual-sensitivity values. Default is \code{7}}
  \item{xminmax}{minimum and maximum limits on the x axis. Default is \code{29}}
  \item{yminmax}{minimum and maximum limits on the y axis. Default is \code{29}}
}

\seealso{\code{\link{saplocmap}}}

\examples{
  stimLoc( "sap", "p24d2", "OD" )
}

\author{Ivan Marin-Franch}