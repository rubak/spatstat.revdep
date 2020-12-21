\name{vfsettings}
\alias{vfsettings}
\docType{data}
\title{Settings of visualField object}
\description{Specifies the structure of a vf object (see \code{\link{vf}}) and information pertaining test patterns and their statistical analysis.}
\usage{data( vfsettings )}
\format{
This structure contains a variable specifying the number of columns with patient and subject data, \code{locini} and several sub-structures with relevant information pertaining different pattern of locations and their analysis. Information exists for the test patterns \code{p24d2}, \code{p30d2}, \code{p10d2}, and \code{sgrnfl}. Each sub-structure contains three items:
\describe{
\item{\code{bs}}{locations that correspond to the ananomical region where the blind spot is. For 24-2, those are locations 26 and 35. There are non for the 10-2}
\item{\code{locnum}}{total number of locations of the testing pattern in which stimuli are presented. For 24-2 there are 54 locations, for 30-2 there are 76, for 10-2 there are 68}
\item{\code{locrPD}}{specifies the rank TD value to be used for the derivation of TD. For 24-2, that ranked location would be 7, corresponding approximately (but not very) to the 85th percentile. For 30-2 the ranked locatoion taken is 10, but this needs fixing???. The way PD is calculated really for 30-2 is by taking the 24-2 locations and finding the seventh largest}
  }
}

\seealso{\code{\link{vf}}}

\keyword{datasets}

\author{Ivan Marin-Franch}