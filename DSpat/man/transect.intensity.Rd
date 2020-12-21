\name{transect.intensity}
\alias{transect.intensity}
\title{Compute expected and observed counts by distance within transect}
\description{
Computes the expected and observed counts for equally-spaced bins of perpendincular distance within
each transect. Expected and observed counts are each a matrix with a row for each transect and a column for each
distance bin.
}
\usage{
transect.intensity(x, epsvu=NULL, obs.ppp, covariates, nclass=NULL, width)
}
\arguments{
  \item{x}{dspat object}
  \item{epsvu}{epsvu setting for fitted model; only uses epsvu[2] value for u}
  \item{obs.ppp}{observation point process}
  \item{covariates}{dataframe of covariates at quadrature points}
  \item{nclass}{number of equally-spaced distance intervals within 0-width/2}
  \item{width}{maximum full transect width over all transects}

}
\value{
  \item{exp.counts}{matrix of expected counts in each distance bin (columns)
                    for each transect (row)}
  \item{obs.counts}{matrix of observed counts in each distance bin (columns)
                    for each transect (row)}
}
\details{
The actual number of distance bins will only match \code{nclass} if it is selected
such that \code{nclass*epsvu[2]} is an even multiple of \code{width/2}.  The
function \code{\link{dspat}} adjusts \code{epsvu[2]} such that it is an even
multiple of \code{width/2} and this function assumes that condition holds.
Sometimes your choice of \code{epsvu[2]} will provide less than optimal choices
for \code{nclass} and in some cases it can only choose a single bin. In these cases,
select another value of \code{epsvu[2]} which is a multiple for \code{width/2}.

If \code{nclass} is not specified then it uses the default of \code{ceiling(sqrt(n))}
intervals.

}
\author{Jeff Laake }

