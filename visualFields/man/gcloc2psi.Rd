\name{gcloc2psi}
\alias{gcloc2psi}
\title{Angle of incidence in the optic nerve head for vf locations}
\description{This is the inverse of the model: it obtains the angle of incidence psi of the average path that passes through position \code{(x,y)}}
\usage{gcloc2psi( xy, r0 = 4 )}
\arguments{
  \item{xy}{Visual field location in Cartesian coordinates}
  \item{r0}{Radius of the optic head nerve. It is a necessity of the model and changing it, changes the calculated average bundle paths. Default value is \code{4}}
}

\value{The angle of incidence psi of the average path that passes through position \code{(x,y)}}

\references{
[1] N. M. Jansonius, J. Nevalainen, B. Selig, L. M. Zangwill, P. A. Sample, W. M. Budde, J. B. Jonas, W. A. Lagreze, P. J. Airaksinen, R. Vonthein, L. A. Levin, J. Paetzold, and U. Schiefer. A mathematical description of nerve fiber bundle trajectories and their variability in the human retina. Vision Research, 49(17):2157-2163, 2009.

[2] N. M. Jansonius, J. Nevalainen, B. Selig, L. M. Zangwill, P. A. Sample, W. M. Budde, J. B. Jonas, W. A. Lagreze, P. J. Airaksinen, R. Vonthein, L. A. Levin, J. Paetzold, and U. Schiefer. Erratum to "a mathematical description of nerve fiber bundle trajectories and their variability in the human retina" [vision research 49(17) (2009) 2157-2163]. Vision Research, 50:1501, 2010.

[3] N. M. Jansonius, J. Schiefer, J. Nevalainen, J. Paetzold, and U. Schiefer. A mathematical model for describing the retinal nerve fiber bundle trajectories in the human eye: Average course, variability, and influence of refraction, optic disc size and optic disc position. Experimental Eye Research, 105:70-78, 2012.
}

\seealso{\code{\link{fiberpathpsi}}}

\examples{
  gcloc2psi( data.frame( x = c( 0, 10 ), y = c( 0, 10 ) ) )
}

\author{Ivan Marin-Franch}