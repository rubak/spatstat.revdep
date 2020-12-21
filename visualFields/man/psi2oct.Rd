\name{psi2oct}
\alias{psi2oct}
\title{Angle of incidence in the OCT scan corresponding to angle of incidence in optic nerve head}
\description{Obtains the angle of incidence in the circular oct scan from the average path that starts at an angle psi from the optic nerve head.}
\usage{psi2oct( psi, diam = 12 )}
\arguments{
  \item{psi}{Angle of incidence of the average bundle path on the optic nerve head}
  \item{diam}{Diameter of the scan. Default value is \code{12} in visual angle, which is the common diameter used in this types of OCT scans}
}

\details{The operation \code{psi2oct( gloc2psi )} maps vf locations with RNFL angles so that the thickness that corresponds to each vf location can be estimated. This can be used in conjuction with \code{vf2gcloc} to correct for ganglion cell body displacement from vf location}

\value{The angle of incidence in the circular oct scan}

\examples{
  psi2oct( 90 + 45 )
}

\author{Ivan Marin-Franch}




