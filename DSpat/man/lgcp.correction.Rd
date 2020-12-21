\name{lgcp.correction}
\alias{lgcp.correction}
\title{Calculate Overdispersion factor for IPP fit via Monte Carlo Integration}
\description{ Calculate Overdispersion factor for IPP fit via Monte Carlo Integration
}
\usage{
lgcp.correction(fit.ppm, fit.lgcp, reps = 100, J.inv, silent = FALSE, lines.psp)
}
\arguments{
  \item{fit.ppm}{fitted model from ppm of spatstat}
  \item{fit.lgcp}{fitted model from lgcp.estK}
  \item{reps}{number of replicates for approximation}
  \item{J.inv}{variance-covariance matrix from fitted ppm model}
  \item{silent}{if FALSE, shows counter for replicates}
  \item{lines.psp}{line segment process}
}
\value{
  \item{J.inv.corr}{Adjusted var-cov matrix}
  \item{u}{score matrix}
}
\author{Devin Johnson}
\seealso{\code{\link{integrate.intensity}}}
