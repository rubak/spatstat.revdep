\name{pdpmap}
\alias{pdpmap}
\title{probability map for pattern deviation }
\description{gets probability map for pattern-deviation values}
\usage{pdpmap( pd )}
\arguments{
  \item{pd}{ pattern-deviation values in dBs }
}
\details{calculates the probability maps for pattern-deviation values}

\value{returns the subjects's pattern-deviation probability map from the age-matched normative reference.}

\seealso{ \code{\link{pdval}}, \code{\link{tdval}}, \code{\link{tdpmap}} }

\examples{
  td  <- tdval( vf91016right )
  pd  <- pdval( td )
  pdp <- pdpmap( pd )
}

\author{Ivan Marin-Franch}