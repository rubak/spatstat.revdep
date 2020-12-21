\name{pdval}
\alias{pdval}
\title{pattern deviation}
\description{gets pattern-deviation values}
\usage{pdval( td )}
\arguments{
  \item{td}{total-deviation values}
}
\details{calculates pattern-deviation values using the (around) 85-th percentile. All in dBs.}

\seealso{ \code{\link{pdpmap}}, \code{\link{tdval}}, \code{\link{tdpmap}} }

\examples{
  td  <- tdval( vf91016right )
  pd  <- pdval( td )
}

\author{Ivan Marin-Franch}