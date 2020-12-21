\name{sdnvghr}
\alias{sdnvghr}
\title{standard deviations of normative values of PD from general-height-rank}
\description{get the standard deviations of PD values from general-height-rank estimates at each location for normative values}
\usage{sdnvghr( vf, smooth = TRUE, smoothFunction = quad2Dfit )}
\arguments{
  \item{vf}{\code{\link{vf}}-object with sensitivity thresholds}
  \item{smooth}{whether to use a function to smooth the results or not. Default is \code{TRUE}}
  \item{smoothFunction}{If \code{smooth} is true, the function to use for smoothing. Default is \code{\link{quad2Dfit}}, a 2D quadratic fit to resulting data. This function is not really a smoothing procedure, but a parametric fit}
}
\details{calculates total-deviation and pattern-deviation values and, from them, their standard deviation at each location}

\value{returns a 3-column data frame with SDs per location. Each row is a location. Column 1 is for sensitivities, column 2 for total deviation, and column 3 for pattern deiation.}

\seealso{ \code{\link{ageLinearModel}}, \code{\link{locperc}}, \code{\link{quad2Dfit}} }

\author{Ivan Marin-Franch}