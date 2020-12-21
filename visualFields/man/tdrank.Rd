\name{tdrank}
\alias{tdrank}
\title{total-deviation rank curve}
\description{get total-deviation rank curve}
\usage{tdrank( td )}
\arguments{
  \item{td}{total deviation values}
}
\details{calculates total-deviation rank curve.}

\seealso{ \code{\link{tdrankperc}} }

\examples{
  td <- tdval( vf91016right )
  tdr <- tdrank( td  )
}

\author{Ivan Marin-Franch}