\name{ghpostd}
\alias{ghpostd}
\title{general height estimated as the rankRef-th ranked TD location}
\description{
gets the general height estimated as the rankRef-th ranked TD location
}
\usage{ghpostd( td, correction = FALSE )}
\arguments{
  \item{td}{vf-object with total-deviation values}
  \item{correction}{obtain general height as difference from mean normal "85th" percentile TD value? Default is \code{FALSE}}
}
\details{calculates the general height estimated as the pos-th ranked TD location}

\value{returns the subject's estimated general height}

\seealso{\code{\link{pdval}}}

\examples{
  td <- tdval( vf91016right[15,] )
  ghseventh <- ghpostd( td )
}

\author{Ivan Marin-Franch}