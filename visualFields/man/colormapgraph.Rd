\name{colormapgraph}
\alias{colormapgraph}
\title{color legend for p-values in td and pd probability plots}
\description{It geneartes color legend for p-values in td and pd probability plots}
\usage{colormapgraph( ncol = 3, mapval = NULL, notSeenAsBlack = TRUE,
                      txtfont = "sans", pointsize = 10,
                      symbol = "circle", size = 1, inch = 0.18  )}
\arguments{
  \item{ncol}{number of columns in where to show the color symbols. Default is \code{3}}
  \item{notSeenAsBlack}{Add a black color-code representing non-seen to the color scheme of the probability maps for td and pd}
  \item{mapval}{map and cutoff values to be used for the generation of the color map. If \code{NULL}, then go to current \code{nv$pmapsettings}. Default is \code{NULL}}
  \item{txtfont}{font of the text with visual-sensitivity values. Default is \code{sans}}
  \item{pointsize}{size of the text with visual-sensitivity values. Default is \code{10}}
  \item{symbol}{The  outer symbol at all locations. Can be any of circles, squares, rectangles, stars. Default is \code{circle}}
  \item{size}{Size of the outer symbol. Default is \code{1}}
  \item{inch}{Maximum size of the outer symbol in inches. Default is \code{0.35}}
}
\details{
it generates a graph with the color legend for p-values in td and pd probability plots given the actual normative values version set for \code{\link{visualFields}}
}

\seealso{\code{\link{vfcolormap}}}

\examples{
  colormapgraph()
}

\author{Ivan Marin-Franch}