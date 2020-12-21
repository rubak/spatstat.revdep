\name{vfplot_sparklines}
\alias{vfplot_sparklines}
\title{Sparklines}
\description{Sparklines to print on top of vfplots}
\usage{vfplot_sparklines( vf, ylim = c( -5, 35 ), collin = NULL, ... )}
\arguments{
  \item{vf}{\code{\link{vf}}-object with sensitivity thresholds, or TD or PD values}
  \item{ylim}{y-limit for the plots for the sparkline.}
  \item{collin}{color for each sparkline to be plotted. Default is \code{NULL}, then all plotted in dark gray}
  \item{...}{Arguments to be passed to methods, such as graphical parameters.}
}
\details{Plots sparklines on top of other vfplots. It is used for the PLR plot in the PoPLR layout.}

\examples{
  \donttest{
    vflayout_poplr( vf91016right )
  }
}

\author{Paul H Artes}