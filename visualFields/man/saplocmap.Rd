\name{saplocmap}
\alias{saplocmap}
\docType{data}
\title{xy-position mapping between HFA device convention and visualFields convention for patterns of location}
\description{
A table with relevant information about test location data for each pattern of locations, 24-2, 10-2, and 30-2, and the conventional Goldman size III stimulus. The G1 pattern is also included for the Goldman size III, size V and size VI stimuli. The convention for visualFields is to use always a right-eye format. That is, a left eye would be "flipped" left-right and location number are counted row-wise from top-left to bottom-right. Information about the size of the stimulus and the corresponding angle of incidence and slope with Jansonious map [1] are included.
}
\usage{data( saplocmap )}
\format{
The structure saplocmap has 1 table for the test pattern \code{p24d2}. The table has six columns:
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

\seealso{\code{\link{fdplocmap}}}

\references{
[1] N. M. Jansonius, J. Nevalainen, B. Selig, L. M. Zangwill, P. A. Sample, W. M. Budde, J. B. Jonas, W. A. Lagreze, P. J. Airaksinen, R. Vonthein, L. A. Levin, J. Paetzold, and U. Schiefer. \emph{A mathematical description of nerve fiber bundle trajectories and their variability in the human retina. Vision Research}, 49, 2009.
}

\author{Ivan Marin-Franch}