\name{vfplot_plr}
\alias{vfplot_plr}
\title{plot with the PLR analysis}
\description{plots the slope values and the corresponding probabiliby category for each location}
\usage{vfplot_plr( sl, pval, vfinfo, newWindow = FALSE,
                   xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL,
                   colorMapType = "pval", colorScale = NULL,
                   txtfont = "sans", pointsize = 10, width = 6,
                   showaxis = FALSE, colaxis = "black" )}
\arguments{
  \item{sl}{slopes}
  \item{pval}{pvalues calculated for the slope}
  \item{vfinfo}{information about the visual field}
  \item{newWindow}{boolean value- TRUE will generate a new window for the plot. Default value is TRUE}
  \item{xmin, xmax, ymin, ymax}{minimum and maximum limits on the x and y axes, in degrees of visual angle. If \code{NULL}, then the limits are the maximum and minimum location values + 2.5 percent of the range. Default is \code{NULL}.}
    \item{colorMapType}{what does color map categorizes. It can be \code{pvals}, \code{slopes}, or years \code{blind}. Default is \code{pvals}}
  \item{colorScale}{Color mapping to use. Default is \code{NULL}. A different default is given depending on \code{colorMapType}}
  \item{txtfont}{font of the text with visual-sensitivity values. Default is \code{sans}}
  \item{pointsize}{size of the text with visual-sensitivity values. Default is \code{10}}
  \item{width}{width of the window. height is calucated using the width, xminmax and yminmax}
  \item{showaxis}{Whether to show axis or not. Default is \code{FALSE}}
  \item{colaxis}{Color of the axis to show, if \code{showaxis} is \code{TRUE}. Default is \code{black}}
}
\details{The vfplot function operates on a single row of visual fileds (vf). Depending on the plot type it will generate a plot representing the sensitivity at each location. The  color scheme elucidates the degree of sensitivity at that location. Locations with zero sensitivity are shown in black. Blind spots have been excluded fronm the plot}

\seealso{\code{\link{poplr}}, \code{\link{poplr_cstat}}, \code{\link{poplr_pstat}}, \code{\link{vflayout_poplr}}}

\references{
[1] N. O'Leary, B. C. Chauhan, and P. H. Artes. \emph{Visual field progression in glaucoma: estimating the overall significance of deterioration with permutation analyses of pointwise linear regression (PoPLR)}. Investigative Ophthalmology and Visual Science, 53, 2012
}

\examples{
  res <- poplr( vf91016right )
  vfplot_plr( res$sl, res$pval, res$vfdata )
}

\author{Ivan Marin-Franch, Paul H Artes, Chaitanya Khadilkar}