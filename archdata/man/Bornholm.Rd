\encoding{UTF-8}
\name{Bornholm}
\alias{Bornholm}
\docType{data}
\title{
Female Iron Age Graves, Bornholm, Denmark
}
\description{
Data on the occurrence of 39 different types of ornamentation in 77 female graves at Iron age sites in in Bornholm, Denmark.
}
\usage{data("Bornholm")}
\format{
  A data frame with 77 observations on the following 42 variables.
  \describe{
    \item{\code{Number}}{Observation Number}
    \item{\code{Site}}{Site/Bural Number}
    \item{\code{Period}}{Chronological period: \code{1a}, \code{1b}, \code{2a}, \code{2b}, \code{2c}, \code{3a}, and \code{3b}}
    \item{\code{N2c}}{count}
    \item{\code{R3d}}{count}
    \item{\code{N2a}}{count}
    \item{\code{Q3b}}{count}
    \item{\code{R3c}}{count}
    \item{\code{N1}}{count}
    \item{\code{Q3c}}{count}
    \item{\code{O1}}{count}
    \item{\code{O2}}{count}
    \item{\code{N2e}}{count}
    \item{\code{I3}}{count}
    \item{\code{R3b}}{count}
    \item{\code{K1a}}{count}
    \item{\code{Q3a}}{count}
    \item{\code{I2}}{count}
    \item{\code{K1c}}{count}
    \item{\code{K1b}}{count}
    \item{\code{H}}{count}
    \item{\code{Q3d}}{count}
    \item{\code{J1d}}{count}
    \item{\code{S1}}{count}
    \item{\code{D}}{count}
    \item{\code{Q2}}{count}
    \item{\code{S3}}{count}
    \item{\code{P2}}{count}
    \item{\code{P4}}{count}
    \item{\code{G3}}{count}
    \item{\code{E2a}}{count}
    \item{\code{P3}}{count}
    \item{\code{R3a}}{count}
    \item{\code{R1}}{count}
    \item{\code{E2b}}{count}
    \item{\code{G2}}{count}
    \item{\code{I1b}}{count}
    \item{\code{G1}}{count}
    \item{\code{F}}{count}
    \item{\code{P1}}{count}
    \item{\code{I1a}}{count}
    \item{\code{A2e}}{count}
  }
}
\details{
Nielsen used data on 39 different types of ornaments from Ørsnes (1966) to seriate a series of 77 Late Germanic Iron Age (CE 550 - 800) graves from Bornholm, Denmark (1988, Table 4 and Figure 7). Baxter re-analyzed the data to illustrate correspondence analysis (1994: 104-107, Table A6). These data were taken from Nielsen's Table 4 showing her seriation. Baxter's version is scrambled in order to evaluate different seriation methods and does not include the ornament types (illustrated in Nielson's Figure 7). The data include Ørsnes's period and subperiod designations (1966).
}
\source{
Baxter, M. J. 1994. \emph{Exploratory Multivariate Analysis in Archaeology}. Edinburgh University Press. Edinburgh.

Nielsen, K. H. 1988. Correspondence Analysis Applied to Hords and Graves of the Germanic Iron Age. In \emph{Multivariate Archaeology: Numerical Approaches in Scandinavian Archaeology}, edited by Torsten Madsen, pp 37-54. Jutland Archaeological Society Publications XXI. Arahus University Press.

Ørsnes, M. 1966.  \emph{Form og stil i Sydskandinaviens yngre germanske jernalder}. Nationalmuseets skrifter. Arkæologisk-historisk række 11. Copenhagen.
}

\examples{
data(Bornholm)
library(MASS)
Bornholm.ca <- corresp(Bornholm[, 4:42], nf=2)
plot(Bornholm.ca$rscore, pch=substring(Bornholm$Period, 1, 1), cex=.75)
boxplot(Bornholm.ca$rscore[, 1]~Bornholm$Period, main="First CA Axis by Period")
}
\keyword{datasets}
