\name{vfsegmentcoord}
\alias{vfsegmentcoord}
\title{calculate line segments to plot in \code{\link{vfplot}} representing the overall orientation of the retinal nerve fibre layer bundles according to Jansonious map [1]}
\description{This function geneartes coordinates for the line segments to be plotted}
\usage{vfsegmentcoord( lineMap, length = 2.5 )}
\arguments{
  \item{lineMap}{patternMap having x,y and slope values}
  \item{length}{length of the line segment in inches. Default is \code{3.75}}
}
\details{This function geneartes coordinates for the line segments to be plotted. Based on the location (x,y) and slope, this function claculates the coordinate set {(x1,y1),(x2,y2)} which is used to plot the line segments}

\references{
[1] N. M. Jansonius, J. Nevalainen, B. Selig, L. M. Zangwill, P. A. Sample, W. M. Budde, J. B. Jonas, W. A. Lagreze, P. J. Airaksinen, R. Vonthein, L. A. Levin, J. Paetzold, and U. Schiefer. \emph{A mathematical description of nerve fiber bundle trajectories and their variability in the human retina. Vision Research}, 49, 2009.
}

\author{Chaitanya Khadilkar, Ivan Marin-Franch}