\name{fiberpathpsi}
\alias{fiberpathpsi}
\title{Average path of a nerve fiber bundle function generator}
\description{Generates a function that renders the average path of a nerve fiber bundle that exits through the optic nerve head with a particular angle}
\usage{fiberpathpsi( psi0, r0 = 4 )}
\arguments{
  \item{psi0}{Angle of incidence of the average bundle path on the optic nerve head}
  \item{r0}{Radius of the optic head nerve. It is a necessity of the model and changing it, changes the calculated average bundle paths. Default value is \code{4}}
}
\details{The function generated works in polar coordinates. The input is the radial coordinate r and output is the angular coordinate. The path in the cartesian \code{(x,y)} space is easily obtained with the function \code{jpolar2cart}}

\value{
A function that returns angular coordinates for each radial coordinate.
}

\references{
[1] N. M. Jansonius, J. Nevalainen, B. Selig, L. M. Zangwill, P. A. Sample, W. M. Budde, J. B. Jonas, W. A. Lagreze, P. J. Airaksinen, R. Vonthein, L. A. Levin, J. Paetzold, and U. Schiefer. A mathematical description of nerve fiber bundle trajectories and their variability in the human retina. Vision Research, 49(17):2157-2163, 2009.

[2] N. M. Jansonius, J. Nevalainen, B. Selig, L. M. Zangwill, P. A. Sample, W. M. Budde, J. B. Jonas, W. A. Lagreze, P. J. Airaksinen, R. Vonthein, L. A. Levin, J. Paetzold, and U. Schiefer. Erratum to "a mathematical description of nerve fiber bundle trajectories and their variability in the human retina" [vision research 49(17) (2009) 2157-2163]. Vision Research, 50:1501, 2010.

[3] N. M. Jansonius, J. Schiefer, J. Nevalainen, J. Paetzold, and U. Schiefer. A mathematical model for describing the retinal nerve fiber bundle trajectories in the human eye: Average course, variability, and influence of refraction, optic disc size and optic disc position. Experimental Eye Research, 105:70-78, 2012.
}

\seealso{\code{\link{gcloc2psi}}}

\examples{
  fiberpathpsi( 90 + 45 )
}

\author{Ivan Marin-Franch}