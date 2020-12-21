\name{vfplot_legoplot}
\alias{vfplot_legoplot}
\title{legoplot for PoPLR analysis}
\description{legoplot for PoPLR analysis}
\usage{vfplot_legoplot( vals, patternMap, vftiles, vfhull,
                        loccolout, loccolin, radius = 2,
                        xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL,
                        txtfont = "sans", pointsize = 10,
                        showaxis = FALSE, colaxis = "black" )}
\arguments{
  \item{vals}{values to show in the lego-plot, typically differences between baseline and final visual field tests}
  \item{patternMap}{locations where to plot symbols, the TD and PD values}
  \item{vftiles}{Tiles for the tessellation calculated with Voronoi polygons}
  \item{vfhull}{Hull defined for the visual field test.}
  \item{loccolout}{color of the polygon surrinding the circle in the lego piece}
  \item{loccolin}{color of the circle in the inside of the lego piece}
  \item{radius}{Radius of the circle inside the lego piece. Default is \code{2}.}
  \item{xmin, xmax, ymin, ymax}{minimum and maximum limits on the x and y axes, in degrees of visual angle. If \code{NULL}, then the limits are the maximum and minimum location values + 2.5 percent of the range. Default is \code{NULL}.}
  \item{txtfont}{font face of the text with visual-sensitivity values. Default is \code{sans}}
  \item{pointsize}{size of the text with visual-sensitivity values. Default is \code{10}}
  \item{showaxis}{Whether to show axis or not. Default is \code{FALSE}}
  \item{colaxis}{Color of the axis to show, if \code{showaxis} is \code{TRUE}. Default is \code{black}}
}
\details{Shows the legoplots}

\examples{
  \donttest{
    vflayout_poplr( vf91016right )
  }
}

\author{Ivan Marin-Franch, Paul H Artes}