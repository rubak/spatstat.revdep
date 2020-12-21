\name{Nelson}
\alias{Nelson}
\docType{data}
\title{
Prehistoric Ceramics at Pueblo San Cristobal, New Mexico, USA
}
\description{
Ceramic distribution in a midden deposit at Pueblo San Cristobal reported by Nels Nelson in 1916.
}
\usage{data(Nelson)}
\format{
  A data frame with 10 observations on the following 8 variables.
  \describe{
    \item{\code{Depth}}{Depth in feet from 1 to 10 for 1 foot arbitrary excavation levels}
    \item{\code{Corrugated}}{Number of corrugated ware ceramics}
    \item{\code{Biscuit}}{Number of Biscuit ware ceramics}
    \item{\code{Type_I}}{Number of two and three color painted ware ceramics}
    \item{\code{Type_II_Red}}{Number of two color glazed red ware ceramics}
    \item{\code{Type_II_Yellow}}{Number of two color glazed yellow ware ceramics}
    \item{\code{Type_II_Gray}}{Number of two color glazed gray ware ceramics}
    \item{\code{Type_III}}{Number of three color glazed ware ceramics}
  }
}
\details{
Data from a midden deposit at San Cristobal in the American Southwest. It has been used as a classic illustration of the potential for creating a relative chronology using frequency seriation of ceramic artifact types. The site was occupied approximately from CE 1350 to 1680. Ceramic artifact fragment counts are presented for each 1-foot (30 cm) arbitrary level excavated in the midden deposit. When converted to percentages (usually excluding the corrugated ware), the data illustrate a classical "battleship curve" like those described in Ford (1962).
}
\source{
Nelson, N. C. 1916. Chronology of the Tano Ruins, New Mexico. \emph{American Anthropologist} 18(2): 159--180.
}
\references{
Carlson, David L. 2017. \emph{Quantitative Methods in Archaeology Using R}. Cambridge University Press, pp 381-386, 390-393.

Ford, J. A. (1962) A Quantitative Method for Deriving Cultural Chronology. \emph{Pan American Union, Technical Manual No} 1.
}
\examples{
data(Nelson)
# Remove Depth and Corrugated and compute percentages
Nelson.pct <- prop.table(as.matrix(Nelson[,3:7]), 1)*100
# Percentages for each type by level
round(Nelson.pct, 2)
# Battleship plot from plotrix package
library(plotrix)
battleship.plot(Nelson.pct)
}
\keyword{datasets}
