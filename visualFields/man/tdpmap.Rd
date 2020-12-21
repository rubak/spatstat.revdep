\name{ tdpmap }
\alias{ tdpmap }
\title{ probability map for total deviation }
\description{gets probability map for total-deviation values}
\usage{tdpmap( td )}
\arguments{
  \item{td}{total-deviation values in dBs}
}
\details{calculates the probability maps for total-deviation values}

\seealso{ \code{\link{tdval}}, \code{\link{pdval}}, \code{\link{pdpmap}} }

\examples{
  td  <- tdval( vf91016right )
  tdp <- tdpmap( td )
}

\author{Ivan Marin-Franch}