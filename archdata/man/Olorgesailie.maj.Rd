\name{Olorgesailie.maj}
\alias{Olorgesailie.maj}
\docType{data}
\title{
Major stone tool classes, Olorgesailie, Kenya
}
\description{
The data represent the number of specimens in each of 6 major artifact classes recovered from 19 localities at the Lower Paleolithic site of Olorgesailie as described in Isaac (1977).
}
\usage{data(Olorgesailie.maj)}
\format{
  A data frame with 19 observations showing the number of specimens for each of 6 stone artifact types.
  \describe{
    \item{\code{Large.cutting.tools}}{Number of large cutting tools}
    \item{\code{Heavy.duty.tools}}{Number of heavy duty tools}
    \item{\code{Large.scrapers}}{Number of large scrapers}
    \item{\code{Other.large.tools}}{Number of other large tools}
    \item{\code{Small.tools}}{Number of small tools}
    \item{\code{Spheroids}}{Number of spheroids}
  }
}
\details{
The data come from Table E1 in Isaac (1977: 239). The \code{rownames} identify localities in the lower, middle and upper strata to provide relative chronological placement. They are in the same order as the columns in the table: LS1 (BBB), LS2 (BBA), LS3(FB), LS4(FB-HL), LS5(FB-I3), MS1a(DE/89 A-L), MS1b(DE/89 A-I), MS2a(DE/89 B-L), MS2b(DE/89 B-I), MS3(DE/89 C), MS4(H/6), MS5(H/9 A), MS6(H/9 AM), MS7(Mid), MS8(Meng), MS9(LHS), US1(TRTrM10), US2(Hog), US3(MFS). Potts (2011) provides updated information on the site complex.
}
\source{
Isaac, Glynn Ll. 1977. \emph{Olorgesailie: Archeological Studies of a Middle Pleistocene Lake Basin in Kenya}. The University of Chicago Press.
}
\references{
Carlson, David L. 2017. \emph{Quantitative Methods in Archaeology Using R}. Cambridge University Press, pp 115-119, 138-142.

Potts, R. 2011. Olorgesailie--Retrospective and current synthesis. In \emph{Casting the net wide: papers in honor of Glynn Isaac and his approach to human origins research}, edited by J. Sept and D. Pilbeam, pp 1--20. American School of Prehistoric Resarch Monographs in Archaeology and Paleoanthropology.
}
\examples{
data(Olorgesailie.maj)
# Chi square after removing the first two columns and simulating the p
# value since there are a number of very small expected values
chisq.test(Olorgesailie.maj, simulate.p.value=TRUE)
# Compute percentages over the localities
Olor.pct <- prop.table(as.matrix(Olorgesailie.maj), 1)*100
boxplot(Olor.pct)
}
\keyword{datasets}
