\name{csplocmap}
\alias{csplocmap}
\docType{data}
\title{xy-position mapping between CSP-custom-device convention (William H Swanson perimeter) and visualFields convention for patterns of location}
\description{
A table with relevant information about test location data for each pattern of locations.
}
\usage{data( csplocmap )}
\format{
The structure saplocmap has 1 table for the test pattern \code{sgrnfl}. The table has six columns:
\describe{
\item{\code{xod}}{stimulus x position}
\item{\code{yod}}{stimulus y position}
\item{\code{loc}}{sequential location number in the original device}
\item{\code{size}}{size of the stimulus presentation}
\item{\code{jmangle}}{angle of incidence in blind spot from Jansonious map}
\item{\code{jmslope}}{orientation of an average bundle at that position of the visual field as calculated from the Jansonious map}
\item{\code{region}}{region of the visual fields in comparison with ONH sector.}
  }
}

\keyword{datasets}

\seealso{\code{\link{saplocmap}},\code{\link{fdplocmap}}}

\references{
[1] N. M. Jansonius, J. Nevalainen, B. Selig, L. M. Zangwill, P. A. Sample, W. M. Budde, J. B. Jonas, W. A. Lagreze, P. J. Airaksinen, R. Vonthein, L. A. Levin, J. Paetzold, and U. Schiefer. \emph{A mathematical description of nerve fiber bundle trajectories and their variability in the human retina. Vision Research}, 49, 2009.
}

\author{Ivan Marin-Franch}