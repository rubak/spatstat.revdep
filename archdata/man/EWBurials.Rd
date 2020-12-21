\name{EWBurials}
\alias{EWBurials}
\docType{data}
\title{
Ernest Witte Cemetery, Austin, County, Texas, U.S.A.
}
\description{
Sex, age, burial group, location, and burial orientation and direction facing from the Ernest Witte site, a Late Archaic cemetery in Texas (Hall 1981).

}
\usage{data(EWBurials)}
\format{
  A data frame with 49 observations on the following 7 variables.
  \describe{
    \item{\code{Group}}{Cemetery group, a factor with levels \code{1}, \code{2}}
    \item{\code{North}}{North grid location of the burial in meters (excavation grid system)}
    \item{\code{West}}{East grid location of the burial in meters (excavation grid system)}
    \item{\code{Age}}{Age category, a factor with levels \code{Fetus}, \code{Infant}, \code{Child}, \code{Adolescent}, \code{Young Adult}, \code{Adult}, \code{Middle Adult}, \code{Old Adult}}
    \item{\code{Sex}}{a factor with levels \code{Female}, \code{Male}}
    \item{\code{Direction}}{circular data in degrees indicating the direction of the individual measured from the head along the vertebral column}
    \item{\code{Looking}}{circular data in degrees indication the direction the individual is facing}
    \item{\code{Goods}}{Presence or absence of grave goods}
  }
}
\details{
The Ernest Witte site in Austin County, Texas contains four burial groups from different time periods. Group 1 includes 60 interments and that occurred between about 2000 and 1200 BCE.  Group 2 is the largest with 148 interments. The burials in this group were interred between about CE 200 and 500. Groups 3 and 4 include only 10 and 13 interments and date to CE 500 to 1500, but are not included in this data set which was taken from Appendix II (Hall 1981). Two of the variables, \code{direction} and \code{looking}, are circular data and require package \code{circular}. Hall (2010) provides a summary of the site and its significance.
}
\source{
Hall, G. D. 1981. Allen's Creek: A Study in the Cultural Prehistory of the Lower Brazos River Valley. \emph{The University of Texas at Austin. Texas Archeological Survey. Texas. Research Report No.} 61.
}
\references{
Carlson, David L. 2017. \emph{Quantitative Methods in Archaeology Using R}. Cambridge University Press, pp 350-357.

Hall, G. D. 2010. Ernest Witte site. \emph{Handbook of Texas Online} \url{http://www.tshaonline.org/handbook/online/articles/bbe05}. Texas State Historical Association.
}
\examples{
data(EWBurials)
xtabs(~Age+Sex+Group, EWBurials)

library(circular)
plot(EWBurials$Direction)
}
\keyword{datasets}
