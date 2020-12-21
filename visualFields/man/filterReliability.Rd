\name{filterReliability}
\alias{filterReliability}
\title{identify visual fields that did not pass filter for reliability}
\description{identify visual fields that did not pass filter for reliability}
\usage{filterReliability( vf, relCriteria = c( 0.2, 0.2, 0.2 ) )}
\arguments{
  \item{vf}{visual field, should have columns perc of false positives, false negatives, and fixation losses}
  \item{relCriteria}{reliability limits for false positives, false negatives, and fixation losses}
}
\value{
returns a list of visual fields that did not pass the reliability criteria
}

\seealso{ \code{\link{lidLensArtifact}} }

\author{Ivan Marin-Franch}