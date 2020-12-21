\name{PitHouses}
\alias{PitHouses}
\docType{data}
\title{
Late Stone Age and Early Sami Iron Age Pithouses in Arctic Norway
}
\description{
The morphology of 45 Arctic Norway pithouses is described in terms of 6 categorical variables.
}
\usage{data("PitHouses")}
\format{
  A data frame with 45 observations on the following 6 variables.
  \describe{
    \item{\code{Hearths}}{a factor with levels \code{None}, \code{One}, \code{Two}, and \code{Charcoal Conc}}
    \item{\code{Depth}}{a factor with levels \code{Deep} and \code{Shallow}}
    \item{\code{Size}}{a factor with levels \code{Small}, \code{Medium}, and \code{Large}}
    \item{\code{Form}}{a factor with levels \code{Oval} and \code{Rectangular}}
    \item{\code{Orient}}{a factor with levels \code{Parallel Coast} and \code{Gabel Toward Coast}}
    \item{\code{Entrance}}{a factor with levels \code{One Side}, \code{Front and One Side}, and \code{None}}
  }
}
\details{
Data on the morphology of pit houses from Arctic Norway described by Engelstad (1988). The data were scanned from Table A7 in Baxter (1994). The category labels are used rather than the numeric values listed in Table A7. The data represent the Group C pithouses as described in Engelstad (1988) which was more variable than groups A or B. The data were converted into an incidence matrix (Table A8 in Baxter (1994) and Table 3 in Englestad (1988)) and used in a multiple correspondence analysis.
}
\source{
Baxter, M. J. 1994. \emph{Exploratory Multivariate Analysis in Archaeology}. Edinburgh University Press.

Engelstad, E. 1988. Pit Houses in Arctic Norway - An Investigation of Their Typology Using Multiple Correspondence Analysis. In \emph{Multivariate Archaeology}, edited by T. Madsen, pp. 71-84. Aarhus University Press.
}
\references{
Carlson, David L. 2017. \emph{Quantitative Methods in Archaeology Using R}. Cambridge University Press, pp 192-197.
}

\examples{
data(PitHouses)
# Crosstabulation of Hearths with Size
PitHouses.tbl <- xtabs(~Hearths+Size, PitHouses)
PitHouses.tbl
barplot(PitHouses.tbl, ylab="Frequency", main="Arctic Norway Pithouses", beside=TRUE,
    legend.text=TRUE, args.legend=list(title="Hearths"))
barplot(prop.table(PitHouses.tbl, 2)*100, ylim=c(0, 60), main="Arctic Norway Pithouses",
    ylab="Percent", beside=TRUE, legend.text=TRUE, args.legend=list(title="Hearths"))
}
\keyword{datasets}
