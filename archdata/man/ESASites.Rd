\encoding{UTF-8}
\name{ESASites}
\alias{ESASites}
\docType{data}
\title{
Early Stone Age Sites - Norway
}
\description{
 Data on 43 Early Stone Age assemblages in Norway come originally from Bølviken et al (1982).
}
\usage{data("ESASites")}
\format{
  A data frame with 43 observations on the following 16 variables.
  \describe{
    \item{\code{TA}}{Tanged Arrows}
    \item{\code{BA}}{Blade Arrows}
    \item{\code{TOA}}{Transverse and Oblique Arrows}
    \item{\code{AA}}{Atypical Arrows}
    \item{\code{M}}{Microliths}
    \item{\code{FK}}{Flake Knives}
    \item{\code{BK}}{Blade Knives}
    \item{\code{NK}}{Notched Knives}
    \item{\code{CFS}}{Core and Flake Scrapers}
    \item{\code{BS}}{Blade Scrapers}
    \item{\code{DS}}{Disc Scrapers}
    \item{\code{Bu}}{Burins}
    \item{\code{Ax}}{Axes}
    \item{\code{Ch}}{Chisels}
    \item{\code{SAx}}{Slate Axes}
    \item{\code{Pf}}{Perforators}
  }
}
\details{
 Data on 43 Early Stone Age (8000 - 4000 BCE) assemblages in Norway come originally from  Bølviken et al (1982). The data were scanned from Table A5 (Appendix A) in Baxter (1994).
}
\source{
Baxter, M. J. 1994. \emph{Exploratory Multivariate Analysis in Archaeology}. Edinburgh University Press. Edinburgh.

Bølviken, E., E. Helskog, K. Helskog, I. M. Holm-Olsen, L. Solheim, and R. Bertelsen. 1982. Correspondence Analysis: An Alternative to Principal Components. \emph{World Archaeology} 14: 41-60.
}
\references{
Carlson, David L. 2017. \emph{Quantitative Methods in Archaeology Using R}. Cambridge University Press, pp 398-410.
}
\examples{
data(ESASites)
NS <- rowSums(ESASites)
NT <- rowSums(ESASites > 0)
plot(NS, NT, xlab="Number of Artifacts", ylab="Number of Types", main="Early Stone Age Sites",
    las=1)
}
\keyword{datasets}
