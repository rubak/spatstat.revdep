\encoding{UTF-8}
\name{Acheulean}
\alias{Acheulean}
\docType{data}
\title{
Seven African Acheulean Sites
}
\description{
Stone tool assemblage data from a paper by Lewis Binford (1972). The sites include Olorgesailie, Isimila, Kalambo Falls, Lochard, Kariandusi, Broken Hill, and Nsongezi. Data include approximate latitude and longitude for each site as well as the frequency for each of 12 stone tool types.
}
\usage{data(Acheulean)}
\format{
   A data frame with 7 observations showing the site location and the number of specimens for each of 12 stone artifact types. The localities are identified by \code{rownames}.
  \describe{
    \item{\code{Lat}}{Latitude (approximate)}
    \item{\code{Long}}{Longitude (approximate)}
    \item{\code{HA}}{Number of handaxes}
    \item{\code{CL}}{Number of cleavers}
    \item{\code{KN}}{Number of knives}
    \item{\code{FS}}{Number of flake scrapers}
    \item{\code{D}}{Number of discoids}
    \item{\code{CS}}{Number of core scrapers}
    \item{\code{P}}{Number of picks}
    \item{\code{CH}}{Number of choppers}
    \item{\code{SP}}{Number of spheroids}
    \item{\code{OLT}}{Number of other large tools}
    \item{\code{SS}}{Number of small scrapers}
    \item{\code{OST}}{Number of other small tools}
  }
}
\details{
Binford (1972) presents the percentages for 12 tool types at 32 assemblages from 7 sites (including Olorgesailie) which was based on Maxine Kleindienst's analysis of Lower Paleolithic Acheulean sites in Africa (1961 and 1962). The data were also analyzed by Glynn Isaac (1977). To create the Acheulean data set, the percentages in the original publication have been converted back to counts by dividing by 100 and multiplying by the number of tools. The assemblages from each site are summed. The largest assemblage is Kalambo Falls with 1349 artifacts and the smallest is Broken Hill (Kabwe) with 94. The \code{rownames} identify each site and an attribute named \code{Variables} provides variable labels for each column.
}
\source{
Binford, L. R. 1972. Contemporary Model Building: Paradigms and the Current State of Paleolithic Research. In \emph{Models in Archaeology}, edited by D. L. Clarke, pp 109--166. Methuen.
}
\references{
Carlson, David L. 2017. \emph{Quantitative Methods in Archaeology Using R}. Cambridge University Press, pp 304-314.

Isaac, Glynn Ll. 1977. \emph{Olorgesailie: Archeological Studies of a Middle Pleistocene Lake Basin in Kenya}. University of Chicago.

Kleindienst, Maxine R. 1961. Variability within the Late Acheulian assemblage in East Africa. \emph{South African Archaeological Bulletin} 16: 35--52.

Kleindienst, Maxine R. 1962. Components of the East African Acheulian assemblage: an analytic approach. In \emph{Actes du IVe Congrès Panafricain de Préhistoire et de l’Étude du Quaternaire}, ed. C. Mortelmans and J. Nenquin, pp 81-105.
}
\examples{
data(Acheulean)
# Compute percentages for each assemblage
Acheulean.pct <- prop.table(as.matrix(Acheulean[,3:14]), 1)*100
round(Acheulean.pct, 2)
plot(OST~HA, Acheulean.pct)
boxplot(Acheulean.pct)
}
\keyword{datasets}
