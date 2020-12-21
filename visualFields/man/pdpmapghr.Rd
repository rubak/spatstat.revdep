\name{pdpmapghr}
\alias{pdpmapghr}
\title{probability map for pattern deviation from global sensitivity estimate}
\description{gets probability map for pattern-deviation values obtained from estimates of general height rank}
\usage{pdpmapghr( pd )}
\arguments{
  \item{pd}{ pattern-deviation values in dBs }
}
\details{calculates the probability maps for pattern-deviation values obtained from estimates from general height rank}

\value{returns the subjects's pattern-deviation probability map from the age-matched normative reference.}

\seealso{ \code{\link{pdval}}, \code{\link{tdval}}, \code{\link{tdpmap}} }

\examples{
  td  <- tdval( vf91016right )
  pd  <- pdvalghr( td )
  pdp <- pdpmapghr( pd )
}

\author{Ivan Marin-Franch}