\name{vfplotloc}
\alias{vfplotloc}
\title{visual field data plot}
\description{plots location-specific data at each location of the visual field. It is the function called by \code{\link{vfplot}}}
\usage{vfplotloc( vals, patternMap, loccol, vftiles, vfhull,
                  xmin, xmax, ymin, ymax,
                  txtfont = "sans", pointsize = 10,
                  showaxis = FALSE, colaxis = "black" )}
\arguments{
  \item{vals}{sensitivity threshold, TD or PD values, or other location-specific values}
  \item{patternMap}{locations where to plot symbols, the TD and PD values}
  \item{loccol}{color to apply to each location. Typically calculated from the normative values}
  \item{vftiles}{Tiles for the tessellation calculated with Voronoi polygons}
  \item{vfhull}{Hull defined for the visual field test.}
  \item{xmin, xmax, ymin, ymax}{minimum and maximum limits on the x and y axes, in degrees of visual angle.}
  \item{txtfont}{font of the text with visual-sensitivity values. Default is \code{sans}}
  \item{pointsize}{size of the text with visual-sensitivity values. Default is \code{10}}
  \item{showaxis}{Whether to show axis or not. Default is \code{FALSE}}
  \item{colaxis}{Color of the axis to show, if \code{showaxis} is \code{TRUE}. Default is \code{black}}
}
\details{The vfplot function operates on a single row of visual fileds (vf). Depending on the plot type it will generate a plot representing the sensitivity at each location. The  color scheme elucidates the degree of sensitivity at that location. Locations with zero sensitivity are shown in black. Blind spots have been excluded fronm the plot}

\seealso{\code{\link{vfplot}}}

\author{Chaitanya Khadilkar, Ivan Marin-Franch}