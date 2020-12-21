\name{tdval}
\alias{tdval}
\title{total deviation}
\description{gets total-deviation values}
\usage{tdval( vf )}
\arguments{
  \item{vf}{ stimulus sensitivities in dBs }
}
\details{calculates the normative reference sensitivities for healthy subjects of the same age as the subject and obtains differences between the visual field sensitivities and mean normal sensitivities. All in dBs.}

\seealso{ \code{\link{tdpmap}}, \code{\link{pdval}}, \code{\link{pdpmap}} }

\examples{
  td <- tdval( vf91016right )
}

\author{Ivan Marin-Franch}