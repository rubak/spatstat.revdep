\name{vflayoutghr}
\alias{vflayoutghr}
\title{printout with results for the visual field}
\description{creates a printout with results for the visual field. This includes sensitivity threshold values with a gray scale plot, and TD and PD values (from general-height-rank estimates) with color coding for the probability map}
\usage{vflayoutghr( vf, pwidth = 8.27, pheight = 11.69, margin = 0.25, filename = NULL,
                    showaxis = FALSE, colaxis = "black" )}
\arguments{
  \item{vf}{\code{\link{vf}}-object with sensitivity thresholds}
  \item{pwidth}{width of the page (in inches). Default is \code{8.27}}
  \item{pheight}{height of the page (in inches). Default is \code{11.69}}
  \item{margin}{margins of the page (in inches). Default is \code{0.25}}
  \item{filename}{file name to save the printout as pdf. If it is saved to pdf, it won't be displayed in the screen. Default is \code{NULL}.}
  \item{showaxis}{Whether to show axis or not. Default is \code{FALSE}}
  \item{colaxis}{Color of the axis to show, if \code{showaxis} is \code{TRUE}. Default is \code{black}}
}

\seealso{\code{\link{vflayout}}, \code{\link{vfplot}}, \code{\link{vfplotloc}}}

\examples{
  \donttest{
    vflayoutghr( vf91016right[15,] )
  }
}

\author{Ivan Marin-Franch}