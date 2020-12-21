\name{fdplocmap}
\alias{fdplocmap}
\docType{data}
\title{xy-position mapping between FDP device convention and visualFields convention for patterns of location}
\description{A table with relevant information about test location data for the pattern of locations 24-2. The convention for visualFields is to use always a right-eye format. That is, a left eye would be "flipped" left-right and location number are counted row-wise from top-left to bottom-right. Information about the size of the stimulus and the corresponding angle of incidence and slope with Jansonious map [1] are included.}
\usage{data( fdplocmap )}
\format{
The structure fdplocmap has 3 tables, one for each test pattern: \code{p24d2}, \code{p10d2}, \code{p30d2}. Each table has six columns:
\describe{
\item{\code{xod}}{stimulus x position}
\item{\code{yod}}{stimulus y position}
\item{\code{loc}}{sequential location number in the original device}
\item{\code{size}}{size of the stimulus presentation}
\item{\code{jmangle}}{angle of incidence in blind spot from Jansonious map}
\item{\code{jmslope}}{orientation of an average bundle at that position of the visual field as calculated from the Jansonious map}
\item{\code{region}}{region of the visual fields in comparison with ONH sector. Garway-Heath map}
  }
}

\keyword{datasets}

\seealso{\code{\link{saplocmap}}}

\references{
[1] N. M. Jansonius, J. Nevalainen, B. Selig, L. M. Zangwill, P. A. Sample, W. M. Budde, J. B. Jonas, W. A. Lagreze, P. J. Airaksinen, R. Vonthein, L. A. Levin, J. Paetzold, and U. Schiefer. \emph{A mathematical description of nerve fiber bundle trajectories and their variability in the human retina. Vision Research}, 49, 2009.
}

\author{Ivan Marin-Franch}