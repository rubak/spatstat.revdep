\name{Mesolithic}
\alias{Mesolithic}
\docType{data}
\title{
British Mesolithic assemblages
}
\description{
Counts of 5 different stone artifact types from 33 Mesolithic sites in Britain.
}
\usage{data("Mesolithic")}
\format{
  A data frame with 33 observations on the following 5 variables.
  \describe{
    \item{\code{Microliths}}{Number of microliths}
    \item{\code{Scrapers}}{Number of scrapers}
    \item{\code{Burins}}{Number of burins}
    \item{\code{Axes}}{Number of axes}
    \item{\code{Saws}}{Number of saws}
  }
}
\details{
Data on 33 Mesolithic (9000 - 4000 BCE) assemblages are a subset Pitts (1979) extracted from a set published by Mellars (1976). The data were scanned from Table A3 (Appendix A) in Baxter (1994).
}
\source{
Baxter, M. J. 1994. \emph{Exploratory Multivariate Analysis in Archaeology}. Edinburgh University Press. Edinburgh.

Mellars, P. 1976. Settlement Patterns and Industrial Variability in the British Mesolithic. In \emph{Problems in Economic and Social Archaeology}, edited by Sieveking, G de G., I. H. Longworth, and K. E. Wilson, pp 375-99. Duckworth, London.

Pitts,  M. W. 1979. Hides and Antlers: A New Look at the Gatherer-Hunter Site at Star Carr, North Yorkshire, England.  \emph{World Archaeology} 11: 32-44.
}
\examples{
data(Mesolithic)
Mesolithic.pct <- prop.table(as.matrix(Mesolithic), 1)*100
apply(Mesolithic.pct, 2, quantile)
cor(Mesolithic.pct)
}
\keyword{datasets}
