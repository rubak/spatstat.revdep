\encoding{UTF-8}
\name{BarmoseI.pp}
\alias{BarmoseI.pp}
\docType{data}
\title{
Piece plotted artifacts from Barmose I, South Zealand, Denmark
}
\description{
Two dimensional locations of 473 artifacts at the Barmose I Maglemosian site used by Blankholm (1991) to illustrate several spatial analysis methods.
}
\usage{data(BarmoseI.pp)}
\format{
  A data frame with 473 observations on the following 4 variables.
  \describe{
    \item{\code{North}}{North coordinate}
    \item{\code{East}}{East coordinate}
    \item{\code{Class}}{Numeric code used by Blankholm: \code{1}, \code{2}, \code{3}, \code{4}, \code{5}, \code{6}, \code{7}, \code{8}, \code{9}, \code{10}, \code{11}}
    \item{\code{Label}}{Artifact type: \code{Scrapers}, \code{Burins}, \code{Lanceolate Microliths}, \code{Microburins}, \code{Flake Axes}, \code{Core Axes}, \code{Square Knives}, \code{Blade/Flake Knives}, \code{Denticulated/Notched Pieces}, \code{Cores}, \code{Core Platforms}}
  }
}
\details{
Barmose I is an early Maglemosian (7500 - 6000 BCE)  site located in Barmosen in South Zealand, Denmark. The site was excavated in 1967-1971 by Axel Johansson (Johansson 1971 and 1990). Flake counts and grid coordinates were taken from Figure 100 in Blankholm (1991) for \code{BarmoseI.grid}. \code{BarmoseI.pp} includes the locations of 473 artifacts from Appendix C of Blankholm's book (1991). 
}
\source{
Blankholm, Hans Peter. 1991. \emph{Intrasite Spatial Analysis in Theory and Practice}. Aarhus University Press.
}
\references{
Carlson, David L. 2017. \emph{Quantitative Methods in Archaeology Using R}. Cambridge University Press, pp 367-377.

Johansson, Axel. 1971. Barmosegruppen. Præboreale Bopladsfund me Skiveøkser i Sydsjælland. Foreløbig Meddelelse. \emph{Historisk Samfund for Præstø Amt}. Årbog 1968, pp. 101--170.

Johansson, Axel. 1990. \emph{Barmosegruppen. Præboreale Bopladsfund i Sydsjælland}. Årbog. Aarhus University Press.
}
\examples{
data(BarmoseI.pp)
plot(North~East, BarmoseI.pp, asp=1, pch=as.numeric(Class))
legend("bottomleft", levels(BarmoseI.pp$Label), pch=1:11, cex=.75)
}
