\name{pdvalghr}
\alias{pdvalghr}
\title{pattern deviation from general height rank (GHr)}
\description{gets the pattern-deviation from the estimator of general height rank (GHr)}
\usage{pdvalghr( td )}
\arguments{
  \item{td}{total-deviation values}
}
\details{calculates pattern-deviation from the estimator of general height rank (GHr). All in dBs.}

\seealso{ \code{\link{pdval}}, \code{\link{pdpmap}}, \code{\link{pdpmapghr}}, \code{\link{tdval}}, \code{\link{tdpmap}} }

\examples{
  td  <- tdval( vf91016right )
  pd  <- pdvalghr( td )
}

\author{Ivan Marin-Franch}