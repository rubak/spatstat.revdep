\name{readGPS}
\alias{readGPS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{GPSbabel read interface}
\description{
  The function reads a data frame from an attached GPS using the external program gpsbabel. The columns of the data frame need to be identified by hand because different GPS order NMEA data in different ways, and the columns should be converted to the correct classes by hand. Once the specifics of a particular GPS are identified, and ways of cleaning erroneous locations are found, the conversion of the output data frame into a usable one may be automated.
}
\usage{
readGPS(i = "garmin", f = "usb:", type="w", invisible=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{i}{INTYPE: a supported file type, default "garmin"}
  \item{f}{INFILE: the appropriate device interface, default "usb:", on Windows for serial interfaces commonly "com4:" or similar}
  \item{type}{"w" waypoints, or "t" track, or others provided in gpsbabel}
  \item{invisible}{Under Windows, do not open an extra window}
  \item{...}{arguments passed through to \code{read.table}}
}
\details{
The function just wraps: gpsbabel -i INTYPE -f INFILE -o tabsep -F - in system(), and reads the returned character vector of lines into a data frame. On some systems, INFILE may not be readable by ordinary users without extra configuration. The gpsbabel program must be present and on the user's PATH for the function to work. Typically, for a given GPS, the user will have to experiment first to find a set of data-cleaning tricks that work, but from then on they should be repeatable.
}
\value{
  A data frame of waypoint values
}
\references{\url{https://www.gpsbabel.org}}
\author{Patrick Giraudoux and Roger Bivand}

\examples{
\dontrun{
#b1 <- readGPS(f="usb:")
#str(b1)
#b2 <- b1[1:172,]
#wp0 <- b2[,c(2,3,4,8,9,19)]
#str(wp0)
#wp0$long <- wp0$V9
#wp0$lat <- as.numeric(as.character(wp0$V8))
#wp0$id <- as.character(wp0$V2)
#wp0$alt <- as.numeric(substring(as.character(wp0$V19), 1,
# (nchar(as.character(wp0$V19))-1)))
#wp0$time <- as.POSIXct(strptime(paste(as.character(wp0$V3),
# as.character(wp0$V4)), format="\%d-\%b-\%y \%H:\%M:\%S"))
#str(wp0)
#wp1 <- wp0[,-(1:6)]
#str(wp1)
#summary(wp1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

