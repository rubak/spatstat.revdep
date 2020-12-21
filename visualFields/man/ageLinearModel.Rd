\name{ageLinearModel}
\alias{ageLinearModel}
\title{linear model for age effect on visual sensitivity}
\description{calculates the linear model to account for mean age effect on visual sensitivity at each location from a sample of controls subjects}
\usage{ageLinearModel( vf, smooth = TRUE, smoothFunction = quad2Dfit )}
\arguments{
  \item{vf}{a \code{\link{vf}} object. It needs to be in a specific format}
  \item{smooth}{whether to use a function to smooth the results or not. Default is \code{TRUE}}
  \item{smoothFunction}{If \code{smooth} is \code{TRUE}, the function to use for smoothing. Default is \code{\link{quad2Dfit}}, a 2D quadratic fit to resulting data. This function is not really a smoothing procedure, but a parametric fit}
}
\details{
The function obtains for each location a regression line of sensitivity threshold on age. The \code{\link{vf}} object may contain data for many visits of many subjects. The number of visits per subject can be variable and hence, to account for that, weighted least-squares linear regression is performed in which the weights for is the inverse of the number of visits for the subject. By default, the slopes and intercepts are "smoothed" by fitting a 2-dimensional quadratic fit, as in [1]. Other \code{smoothFunctions} can be defined instead using other type of parametric or non-parametric fits. The 2-dimensional quadratic fit has been found to be innapropriate for the stimulus used in frequency-doubling perimetry (FDP) [2].

All the data passed to the function must belong to the same perimeter \code{tperimetry}, the same pattern of locations \code{talgorithm}, and the same presentation algorithm \code{tpattern}.
}

\value{
returns a table with intercepts and slopes modelling age-related mean sensitivity change.
}

\references{
[1] A. Heijl, G. Lindgren, and J. Olsson. \emph{A package for the statistical analysis of visual fields}. Documenta Ophthalmologica Proceedings Series, 49, 1987

[2] A. J. Anderson, C. A. Johnson, M. Fingeret, J. L. Keltner, P. G. D. Spry, M. Wall, and J. S. Werner. \emph{Characteristics of the normative database for the humphrey matrix perimeter}. Investigative Ophthalmology and Visual Science, 46, 2005
}

\seealso{\code{\link{sdnv}}, \code{\link{locperc}}, \code{\link{quad2Dfit}}}

\author{Ivan Marin-Franch}