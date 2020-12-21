\name{lidLensArtifact}
\alias{lidLensArtifact}
\title{visual fields with lid or lens artifacts}
\description{identifies in a very coarse way visual fields that may have been affected by lid or lens artifacts. This function is to be used only with control data. It is probably better not to use it and look to the visual-fields printouts directly to identify those with lid and lens artifacts}
\usage{lidLensArtifact( vf, min_dB = c( 12 ) )}
\arguments{
  \item{vf}{visual-field object}
  \item{min_dB}{lowest dB value considered normal}
}

\value{
it returns a list of indices of visual fields in \code{\link{vf}} that may have been affected by lid or lens artifact
}

\seealso{ \code{\link{filterReliability}} }

\author{Ivan Marin-Franch}