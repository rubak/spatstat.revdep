\name{BACups}
\alias{BACups}
\docType{data}
\title{
Bronze Age Cups from Italy
}
\description{
Measurements on Early and Late Bronze Age ceramic cups from Italy analyzed by Lukesh and Howe (1978).
}
\usage{data("BACups")}
\format{
  A data frame with 60 observations on the following 6 variables.
  \describe{
    \item{\code{RD}}{Rim Diameter}
    \item{\code{ND}}{Neck Diameter}
    \item{\code{SD}}{Shoulder Diameter}
    \item{\code{H}}{Total Height}
    \item{\code{NH}}{Neck Height}
    \item{\code{Phase}}{Chronological Phase: Protoapennine, Subapennine}
  }
}
\details{
These data on the dimensions of Bronze Age cups from Italy are a subset extracted from a set published by Lukesh and Howe  (1978) of the specimens for which full data was available. The data were scanned from Table A4 (Appendix A) in Baxter (1994). The Protoapennine cups are Early Bronze Age while the Subapennine cups are Late Bronze Age.
}
\source{
Baxter, M. J. 1994. \emph{Exploratory Multivariate Analysis in Archaeology}. Edinburgh University Press. Edinburgh.

Lukesh S. S. and S. Howe 1978. Protoapennine vs. Subapennine: Mathematical Distinction Between Two Ceramic Phases. \emph{Journal of Field Archaeology} 5: 339-47.
}
\examples{
data(BACups)
by(BACups[, -6], BACups$Phase, summary)
plot(RD~H, BACups, pch=as.numeric(Phase))
legend("topleft", levels(BACups$Phase), pch=1:2)
}
\keyword{datasets}
