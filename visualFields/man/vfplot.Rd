\name{vfplot}
\alias{vfplot}
\title{visual field plot}
\description{plots location-specific sensitivity thresholds, or TD or PD values at each location of the visual field with the corresponding color code for the probability map}
\usage{vfplot( vf, plotType = "vf",
               xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL,
               notSeenAsBlack = TRUE, newWindow = FALSE,
               txtfont = "sans", pointsize = 10, width = 6,
               showaxis = FALSE, colaxis = "black" )}
\arguments{
  \item{vf}{\code{\link{vf}}-object with sensitivity thresholds, or TD or PD values}
  \item{plotType}{Type of plot to show. It can be visual sensitivities (\code{vf}), total-deviation values (\code{td}), or pattern-deviation values (\code{pd}). Default is \code{vf}}
  \item{xmin, xmax, ymin, ymax}{minimum and maximum limits on the x and y axes, in degrees of visual angle. If \code{NULL}, then the limits are the maximum and minimum location values + 2.5 percent of the range. Default is \code{NULL}.}
  \item{notSeenAsBlack}{whether to plot non-seen as black. Default is \code{TRUE}}
  \item{newWindow}{whether to generate a new window for the plot. This becomes important when working with RStudio. Default value is \code{TRUE}}
  \item{txtfont}{font face of the text with visual-sensitivity values. Default is \code{sans}}
  \item{pointsize}{size of the text with visual-sensitivity values. Default is \code{10}}
  \item{width}{width of the window. height is calculated using the width, xmin, xmax, ymin, and ymax. Default is \code{6}}
  \item{showaxis}{Whether to show axis or not. Default is \code{FALSE}}
  \item{colaxis}{Color of the axis to show, if \code{showaxis} is \code{TRUE}. Default is \code{black}}
}
\details{This function 
The vfplot function operates on a single row of visual fileds (vf). Depending on the plot type it will generate a plot representing the sensitivity at each location. The  color scheme elucidates the degree of sensitivity at that location. Locations with zero sensitivity are shown in black. Blind spots have been excluded fronm the plot}

\seealso{\code{\link{vfplotloc}}}

\examples{
  vfplot( vf91016right[15,], plotType = "td" )
}

\author{Chaitanya Khadilkar, Ivan Marin-Franch}