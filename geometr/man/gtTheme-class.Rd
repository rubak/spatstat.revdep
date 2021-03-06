% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gtTheme.R
\docType{class}
\name{gtTheme-class}
\alias{gtTheme-class}
\alias{themeClass}
\title{Theme class (S4) and methods}
\description{
An \code{gtTheme} stores a theme to \code{\link{visualise}} vector and raster
objects. It is recommended to use \code{\link{setTheme}} to modify a
\code{gtTheme}, because it carries out all the checks and makes sure that
names of the parameters are properly matched.
}
\section{Slots}{

\describe{
\item{\code{title}}{[\code{named list(3)}]\cr properties of the title.}

\item{\code{box}}{[\code{named list(4)}]\cr properties of the bounding box.}

\item{\code{xAxis}}{[\code{named list(5)}]\cr properties of the x-axis, its labels
and ticks.}

\item{\code{yAxis}}{[\code{named list(5)}]\cr properties of the y-axis, its labels
and ticks.}

\item{\code{grid}}{[\code{named list(5)}]\cr properties of the major and minor grid.}

\item{\code{legend}}{[\code{named list(10)}]\cr properties of the legend, its title,
labels, ticks and bounding box.}

\item{\code{vector}}{[\code{named list(7)}]\cr properties of a vector object.}

\item{\code{raster}}{[\code{named list(2)}]\cr properties of a raster object.}
}}

