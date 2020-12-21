\name{integrate.intensity}
\alias{integrate.intensity}
\title{Integrated intensity of fitted model}
\description{
Compute intensity and its integration (abundance) and measures of precision with and
without over-dispersion correction}
\usage{
integrate.intensity(x, dimyx=NULL, eps=NULL, se=FALSE, od=FALSE,
                    reps=100, silent=FALSE, J.inv=NULL, showplot=TRUE)
}
\arguments{
  \item{x}{dspat object}
  \item{dimyx}{number of y,x pixels}
  \item{eps}{height and width of pixels}
  \item{se}{if TRUE, compute std error of abundance and log-normal ci}
  \item{od}{if TRUE and se=TRUE, also compute over-dispersion corrected std error of abundance and log-normal ci}
  \item{reps}{number of reps for MC integration for over-dispersion correction}
  \item{silent}{if FALSE, show progress on MC integration}
  \item{J.inv}{var-cov matrix from fitted model}
  \item{showplot}{if TRUE show Poisson and empirical and fitted K-functions}
}
\value{
  \item{Abundance}{Estimate of expected abundance in the study area}
  \item{distribution}{dataframe containing N (predicted number of points in the cell),x,y (x,y coordinates of cell) and covariates used in the model}
  \item{precision}{List containing se, lcl.95, ucl.95, J.inv, and b.vec}
  \item{precision.od}{For over-dispersion estimate a list containing se, lcl.95, ucl.95, J.inv, and b.vec}
  \item{lambda}{estimated intensity image}
  \item{W}{window mask for study area}
}
\details{
Either \code{dimyx} or \code{eps} can be specified.  If neither specified then
it uses the first covariate image in the dspat object to set the intensity grid.  If more
than one are specified then others are ignored with their priority for use matching
the order they are listed above.
}
\author{Devin Johnson; Jeff Laake}
\seealso{\code{\link{lgcp.correction}}}

