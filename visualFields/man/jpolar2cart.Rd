\name{jpolar2cart}
\alias{jpolar2cart}
\title{converts to \code{(x,y)} in degrees from polar coordinates}
\description{converts to \code{(x,y)} in degrees from polar coordinates. It is the inverse of \code{( cart2jpolar)}}
\usage{jpolar2cart( rpsi )}
\arguments{
  \item{rpsi}{Visual field location in polar coordinates}
}
\details{Input rpsi needs to be a data frame. It returns a data frame with the Cartesian coordinates}

\value{
Data frame with the (x,y) Cartesian coordinates
}

\seealso{\code{\link{cart2jpolar}}}

\examples{
  jpolar2cart( cart2jpolar( data.frame( x = c( 0, 10 ), y = c( 0, 10 ) ) ) )
}

\author{Ivan Marin-Franch}