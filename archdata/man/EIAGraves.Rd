\name{EIAGraves}
\alias{EIAGraves}
\docType{data}
\title{
Early Iron Age Graves - Tell el-Far'ah, Palestine
}
\description{
Counts of 52 different ceramic types in 6 large tombs and 10 broadly contemporaneous groups of tombs.
}
\usage{data("EIAGraves")}
\format{
  A data frame with 52 rows (ceramic types) found in 16 units (a grave or a group of graves).
  \describe{
    \item{\code{Type}}{Ceramic type number}
    \item{\code{G100}}{19 broadly contemporaneous graves and tombs}
    \item{\code{G200B}}{30 broadly contemporaneous graves and tombs}
    \item{\code{G200C}}{28 broadly contemporaneous graves and tombs}
    \item{\code{G201}}{An indidivual tomb}
    \item{\code{G229}}{An indidivual tomb}
    \item{\code{G500N}}{19 broadly contemporaneous graves and tombs}
    \item{\code{G532}}{An indidivual tomb}
    \item{\code{G542}}{An indidivual tomb}
    \item{\code{G552}}{An indidivual tomb}
    \item{\code{G562}}{An indidivual tomb}
    \item{\code{G600}}{52 broadly contemporaneous graves and tombs}
    \item{\code{G800}}{39 broadly contemporaneous graves and tombs}
    \item{\code{G900B}}{41 broadly contemporaneous graves and tombs}
    \item{\code{G900L}}{3 broadly contemporaneous graves and tombs}
    \item{\code{G900S}}{5 broadly contemporaneous graves and tombs}
    \item{\code{G900U}}{7 broadly contemporaneous graves and tombs}
  }
}
\details{
The data on counts of 52 different ceramic types in 6 large tombs and 10 broadly contemporaneous groups of tombs come from Tell el-Far'ah (South), Palestine. They were originally published in McClellan (1979).  The data were scanned from Table 2.5 in Baxter (2003, p. 25-6). The 52 rows correspond to different pottery types found in association with the burials.
}
\source{
Baxter, M. J. 2003. \emph{Statistics in Archaeology}. Arnold, London.

McClellan, T. L. 1979. Chronology of the 'Philistine' Burials at Tell el-Farah (South). \emph{Journal of Field Archaeology} 6: 57-73.
}
\examples{
data(EIAGraves)
# How many ceramics of each type?
# Exclude the first column which is the ceramic type number
rowSums(EIAGraves[, -1])
# How many tomb groups contain each type?
rowSums(EIAGraves[, -1]>0)
# How many ceramics in each tomb group?
colSums(EIAGraves[, -1])
# How many types are found in each tomb group?
colSums(EIAGraves[, -1]>0)
}
\keyword{datasets}
