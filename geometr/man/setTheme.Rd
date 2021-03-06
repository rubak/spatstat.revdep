% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualise.R
\name{setTheme}
\alias{setTheme}
\title{Create a new theme}
\usage{
setTheme(
  from = NULL,
  title = NULL,
  box = NULL,
  xAxis = NULL,
  yAxis = NULL,
  grid = NULL,
  legend = NULL,
  scale = NULL,
  vector = NULL,
  raster = NULL
)
}
\arguments{
\item{from}{[\code{gtTheme}]\cr an gtTheme object.}

\item{title}{[\code{named list(.)}]\cr \code{plot = TRUE/FALSE},
\code{fontsize} and \code{colour} of the title.}

\item{box}{[\code{named list(.)}]\cr \code{plot = TRUE/FALSE},
\code{linewidth}, \code{linetype} and \code{linecol} of the bounding box
(not supported recently).}

\item{xAxis}{[\code{named list(.)}]\cr \code{plot = TRUE/FALSE}, number of
\code{bins} and \code{margin} of the x-axis,\cr\cr label [\code{named
list(.)}]\cr \code{plot = TRUE/FALSE}, \code{title}, \code{fontsize},
\code{colour} and \code{rotation} of the x-axis label,\cr\cr ticks
[\code{named list(.)}]\cr \code{plot = TRUE/FALSE}, \code{fontsize},
\code{colour} and number of \code{digits} to which to round the x-axis
ticks.}

\item{yAxis}{[\code{named list(.)}]\cr \code{plot = TRUE/FALSE}, number of
\code{bins} and \code{margin} of the y-axis,\cr\cr label [\code{named
list(.)}]\cr \code{plot = TRUE/FALSE}, \code{title}, \code{fontsize},
\code{colour} and \code{rotation} of the y-axis label,\cr\cr ticks
[\code{named list(.)}]\cr \code{plot = TRUE/FALSE}, \code{fontsize},
\code{colour} and number of \code{digits} to which to round the y-axis
ticks.}

\item{grid}{[\code{named list(.)}]\cr \code{plot = TRUE/FALSE},
\code{colour}, \code{linetype} and \code{linewidth} of the major and minor
grid and whether or not to plot the \code{minor = TRUE/FALSE} grid.}

\item{legend}{[\code{named list(.)}]\cr \code{plot = TRUE/FALSE}, number of
\code{bins}, \code{ascending = TRUE/FALSE} order of values and the
\code{sizeRatio} of plot and legend, \cr\cr label [\code{named list(.)}]\cr
\code{plot = TRUE/FALSE}, \code{fontsize} and \code{colour} of the legend
labels, \cr\cr box [\code{named list(.)}]\cr \code{plot = TRUE/FALSE},
\code{linetype}, \code{linewidth} and \code{colour} of the legend box.}

\item{scale}{[\code{named list(.)}]\cr \code{param = 'someParameter'} and
\code{to = 'someAttribute'} to which to scale 'someParameter' to. Whether
or not to use the values' \code{identity}, the value \code{range} that
shall be represented by the scale and the number of \code{bins}.}

\item{vector}{[\code{named list(.)}]\cr \code{linecol}, \code{fillcol},
\code{missingcol}, \code{linetype}, \code{linewidth}, \code{pointsize} and
\code{pointsymbol} of a vector object.}

\item{raster}{[ \code{named list(.)}]\cr \code{fillcol} of a raster object.}
}
\description{
Assign parameters in a gtTheme to create a new theme.
}
\examples{
input <- gtRasters$continuous
(myTheme <- setTheme(title = list(plot = FALSE)))

visualise(input, theme = myTheme)
}
