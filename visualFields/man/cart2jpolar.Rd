\name{cart2jpolar}
\alias{cart2jpolar}
\title{convert from Cartesian coordinate to polar coordinates for use with the Jansonius map}
\description{converts from \code{(x,y)} in degrees to polar coordinates. This is a necessary step to compute average fiber paths as a function of their angle of incidence in the optic nerve head with the Jansonius map}
\usage{cart2jpolar( xy )}
\arguments{
  \item{xy}{Visual field location in Cartesian coordinates}
}
\details{Input xy needs to be a data frame. It returns a data frame with the radial and angular coordinates}

\value{
Data frame with the radial and angular coordinates
}

\seealso{\code{\link{jpolar2cart}}}

\examples{
  cart2jpolar( data.frame( x = c( 0, 10 ), y = c( 0, 10 ) ) )
}

\author{Ivan Marin-Franch}