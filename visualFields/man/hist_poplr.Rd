\name{hist_poplr}
\alias{hist_poplr}
\title{histogram for PoPLR analysis}
\description{plots an histogram of Fisher S statistic combining all location p-values}
\usage{hist_poplr( s, pval, sp, nLoc = 52, showLabels = TRUE,
                   txtfont = "sans", pointsize = 10,
                   plt = c( 0, 1.0, 0.3, 1 ), mpg = c( 1.85, 0.5, 0 ),
                   coltxt  = rgb( red = 1.0, green = 0.0, blue = 0.0 ),
                   colhist = rgb( red = 1.0, green = 0.0, blue = 0.0, alpha = 0.5 ))}
\arguments{
  \item{s}{observed S-statistic for the pointwise significance tests that slopes are lower than (or greater than) a specified reference (typically zero)}
  \item{pval}{p-value for the pointwise significance tests that slopes are lower than (or greater than) a specified reference (typically zero)}
  \item{sp}{S-statistics calculated for the permutation of visual fields}
  \item{nLoc}{number of locations in the visual field to be analyzed. For 24-2, it is 52 (54 minus the locations in the blind spot). Default is \code{52}.}
  \item{showLabels}{Whether or not to show histogram x-labels. Default is \code{TRUE}.}
  \item{plt}{Graphics parameter plt. Default is \code{c( 0, 1.0, 0.3, 1 )}}
  \item{mpg}{Graphics parameter mpg. Default is \code{c( 1.85, 0.5, 0 )}}
  \item{txtfont}{font of the text with visual-sensitivity values. Default is \code{sans}}
  \item{pointsize}{size of the text with visual-sensitivity values. Default is \code{10}}
  \item{colhist}{Histogram color. Default is \code{rgb( red = 1.0, green = 0.0, blue = 0.0 )}}
  \item{coltxt}{Text color. Default is \code{rgb( red = 1.0, green = 0.0, blue = 0.0 )}}
}

\seealso{\code{\link{poplr}}, \code{\link{poplr_cstat}}, \code{\link{poplr_pstat}}, \code{\link{vflayout_poplr}}, \code{\link{vfplot_plr}}}

\references{
[1] N. O'Leary, B. C. Chauhan, and P. H. Artes. \emph{Visual field progression in glaucoma: estimating the overall significance of deterioration with permutation analyses of pointwise linear regression (PoPLR)}. Investigative Ophthalmology and Visual Science, 53, 2012
}

\author{Ivan Marin-Franch}