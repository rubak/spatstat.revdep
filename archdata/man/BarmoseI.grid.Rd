\encoding{UTF-8}
\name{BarmoseI.grid}
\alias{BarmoseI.grid}
\docType{data}
\title{
Flakes per grid unit from Barmose I, South Zealand, Denmark
}
\description{
Flake counts for each of 107 contiguous grid units at the Barmose I Maglemosian site used by Blankholm (1991) to illustrate several spatial analysis methods.
}
\usage{data(BarmoseI.grid)}
\format{
  A data frame with 107 observations on the following 3 variables.
  \describe{
    \item{\code{North}}{North coordinate of southwest corner of unit}
    \item{\code{East}}{East coordinate of southwest corner of unit}
    \item{\code{Debitage}}{Number of flakes}
  }
}
\details{
Barmose I is an early Maglemosian (7500 - 6000 BCE)  site located in Barmosen in South Zealand, Denmark. The site was excavated in 1967-1971 by Axel Johansson (Johansson 1971 and 1990). Flake counts and grid coordinates were taken from Figure 100 in Blankholm (1991) for \code{BarmoseI.grid}. \code{BarmoseI.pp} includes the locations of 473 artifacts from Appendix C of Blankholm's book.
}
\source{
Blankholm, Hans Peter. 1991. \emph{Intrasite Spatial Analysis in Theory and Practice}. Aarhus University Press.
}
\references{
Carlson, David L. 2017. \emph{Quantitative Methods in Archaeology Using R}. Cambridge University Press, pp 358-367.

Johansson, Axel. 1971. Barmosegruppen. Præboreale Bopladsfund me Skiveøkser i Sydsjælland. Foreløbig Meddelelse. \emph{Historisk Samfund for Præstø Amt}. Årbog 1968, pp. 101-170.

Johansson, Axel. 1990. \emph{Barmosegruppen. Præboreale Bopladsfund i Sydsjælland}. Årbog. Aarhus University Press.
}
\examples{
data(BarmoseI.grid)
plot(North~East, BarmoseI.grid, xlim=c(0, 12), ylim=c(0, 14), type="n", asp=1)
with(BarmoseI.grid, text(East+.5, North+.5, Debitage, cex=.8))
}
\keyword{datasets}
