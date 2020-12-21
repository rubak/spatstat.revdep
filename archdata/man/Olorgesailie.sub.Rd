\encoding{UTF-8}
\name{Olorgesailie.sub}
\alias{Olorgesailie.sub}
\docType{data}
\title{
Stone tool subclasses, Olorgesailie, Kenya
}
\description{
The data represent the number of specimens in each of 16 artifact subclasses recovered from 19 localities at the Lower Paleolithic site of Olorgesailie as described in Isaac (1977).
}
\usage{data(Olorgesailie.sub)}
\format{
  A data frame with 19 observations showing the stratum, locality and the number of specimens for each of 16 stone artifact types.
  \describe{
    \item{\code{Strat}}{stratum: \code{Lower}, \code{Middle}, \code{Upper}}
    \item{\code{Locality}}{Locality}
    \item{\code{HA}}{Number of handaxes}
    \item{\code{PHA}}{Number of pick-like handaxes}
    \item{\code{CHA}}{Number of chisel handaxes}
    \item{\code{CL}}{Number of cleavers}
    \item{\code{KN}}{Number of knives}
    \item{\code{BLCT}}{Number of broken large cutting tools}
    \item{\code{PAT}}{Number of picks and trièdres}
    \item{\code{CH}}{Number of choppers}
    \item{\code{CS}}{Number of core scrapers}
    \item{\code{LFS}}{Number of large flake scrapers}
    \item{\code{CB}}{Number of core bifaces}
    \item{\code{OLT}}{Number of other large tools}
    \item{\code{SSS}}{Number of small scrapers simple}
    \item{\code{SSNP}}{Number of small scrapers nosed point}
    \item{\code{OST}}{Number of other small tools}
    \item{\code{SP}}{Number of spheroids}
  }
}
\details{
The data come from Table E1 in Isaac (1977: 239). The \code{Locality} contains the column headings in the original table. The \code{rownames} are the same as those in \code{Olorgesailie.maj}. The attribute \code{Variables} in the data frame includes the full variable names. Potts (2011) provides updated information on the site complex.
}
\source{
Isaac, Glynn Ll. 1977. \emph{Olorgesailie: Archeological Studies of a Middle Pleistocene Lake Basin in Kenya}. The University of Chicago Press.
}
\references{
Carlson, David L. 2017. \emph{Quantitative Methods in Archaeology Using R}. Cambridge University Press, pp 280-293.

Potts, R. 2011. Olorgesailie--Retrospective and current synthesis. In \emph{Casting the net wide: papers in honor of Glynn Isaac and his approach to human origins research}, edited by J. Sept and D. Pilbeam, pp 1--20. American School of Prehistoric Research Monographs in Archaeology and Paleoanthropology.
}
\examples{
data(Olorgesailie.sub)
# Chi square after removing the first two columns and simulating the p
# value since there are a number of very small expected values
chisq.test(Olorgesailie.sub[,3:18], simulate.p.value=TRUE)
# Compute percentages over the localities
Olor.pct <- prop.table(as.matrix(Olorgesailie.sub[,3:18]), 1)*100
boxplot(Olor.pct, cex.axis=.7)
}
\keyword{datasets}
