\name{vflayout_progress}
\alias{vflayout_progress}
\title{progression analysis layout}
\description{A layout showing graphs and stats for progression analysis.}
\usage{vflayout_progress( vf, plotType, grp = 2, nperm = 5000,
                          colorMapType = "pval", colorScale = NULL,
                          filename = NULL,
                          pwidth = 8.27, pheight = 11.69, margin = 0.25,
                          showaxis = FALSE, colaxis = "black" )}
\arguments{
  \item{vf}{\code{\link{vf}} object with threshold sensitivities, td, or pd values}
  \item{plotType}{Type of plot to show. It can be visual sensitivities (\code{vf}), total-deviation values (\code{td}), or pattern-deviation values (\code{pd})}
  \item{grp}{how many visual fields to group. Default is \code{3}}
  \item{nperm}{number of permutations. Default is \code{5000}}
  \item{colorMapType}{what does color map categorizes. It can be \code{pvals}, \code{slopes}, or years \code{blind}. Default is \code{pvals}}
  \item{colorScale}{Color mapping to use. Default is \code{NULL}. A different default is given depending on \code{colorMapType}}
  \item{filename}{file name to save the printout as pdf. If it is saved to pdf, it won't be displayed in the screen. Default is \code{NULL}.}
  \item{pwidth}{width of the page (in inches). Default is \code{8.27}}
  \item{pheight}{height of the page (in inches). Default is \code{11.69}}
  \item{margin}{margins of the page (in inches). Default is \code{0.25}}
  \item{showaxis}{Whether to show axis or not. Default is \code{FALSE}}
  \item{colaxis}{Color of the axis to show, if \code{showaxis} is \code{TRUE}. Default is \code{black}}
}

\seealso{\code{\link{vflayout}}, \code{\link{vflayout_poplr}}}

\examples{
  \donttest{
    vflayout_progress( vf91016left, plotType = "td" )
  }
}

\author{Ivan Marin-Franch}