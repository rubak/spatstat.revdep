\name{vfcolormap}
\alias{vfcolormap}
\title{Color code for TD or PD probability maps}
\description{returns the RGB values specifying the color code for the TD or PD probability maps}
\usage{vfcolormap( map, mapval = NULL )}
\arguments{
  \item{map}{array with the percentile for TD or PD probability map}
  \item{mapval}{map and cutoff values to be used for the generation of the color map. If \code{NULL}, then go to current \code{nv$pmapsettings}. Default is \code{NULL}}
}

\seealso{\code{\link{vfplot}}, \code{\link{vfgrayscale}}}

\author{Chaitanya Khadilkar, Ivan Marin-Franch}