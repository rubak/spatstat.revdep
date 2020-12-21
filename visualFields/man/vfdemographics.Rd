\name{vfdemographics}
\alias{vfdemographics}
\title{demographics and statistics of sample in visual-fields object vf}
\description{gets demographics and weighted statistics of sample in visual-fields object vf}
\usage{vfdemographics( vf )}
\arguments{
  \item{vf}{A \code{\link{vf}}-object with sensitivity thresholds, total-deviation or pattern-deviation values}
}

\value{returns a table with demographics and weighted statistics}

\examples{
  vfdemog <- vfdemographics( vf91016right )
}

\author{Ivan Marin-Franch}