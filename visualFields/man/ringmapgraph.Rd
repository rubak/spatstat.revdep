\name{ringmapgraph}
\alias{ringmapgraph}
\title{ring legend for p-values}
\description{It geneartes ring legend for p-values}
\usage{ringmapgraph( ncol = 3, mapval = NULL, txtfont = "mono", pointsize = 7,
               outerSymbol = "circles", innerSymbol = "circles",
               outerSize = 1, innerSize = 1,
               outerInch = 0.2, innerInch = 0.1,
               outerBorderThickness = 2, innerBorderThickness = 2 ) }
\arguments{
  \item{ncol}{number of columns in where to show the color symbols. Default is \code{3}}
  \item{mapval}{map and cutoff values to be used for the generation of the color map. Default is \code{NULL}}
  \item{txtfont}{font of the text with visual-sensitivity values. Default is \code{mono}}
  \item{pointsize}{size of the text with visual-sensitivity values. Default is \code{7}}
  \item{outerSymbol}{The  outer symbol at all locations. Can be any of circles, squares, rectangles, stars. Default is \code{circle}}
  \item{innerSymbol}{The  inner symbol at all locations. Can be any of circles, squares, rectangles, stars. Default is \code{circle}}
  \item{outerSize}{Size of the outer symbol. Default is \code{1}}
  \item{innerSize}{Size of the inner symbol. Default is \code{1}}
  \item{outerInch}{Maximum size of the outer symbol in inches. Default is \code{0.2}}
  \item{innerInch}{Maximum size of the inner symbol in inches. Default is \code{0.1}}
  \item{outerBorderThickness}{Thickness of outer border for outer symbol denoting statistical significance. Default is \code{2}}
  \item{innerBorderThickness }{Thickness of the inner border that represents statistical significance. Default is \code{2}}
}

\seealso{\code{\link{vfcolormap}}}

\author{Ivan Marin-Franch}