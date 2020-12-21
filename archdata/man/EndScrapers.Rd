\name{EndScrapers}
\alias{EndScrapers}
\docType{data}
\title{
Upper Paleolithic End Scrapers from Castenet A and Ferrassie H, France
}
\description{
Data on 3000 Upper Paleolithic end scrapers from two sites analyzed by James Sackett (1966) and reanalyzed by Dwight Read (1974 and 2007). 
}
\usage{data(EndScrapers)}
\format{
  A data frame with 48 observations on the following 6 variables.
  \describe{
    \item{\code{Width}}{Width: \code{Narrow}, \code{Wide}}
    \item{\code{Sides}}{Sides: \code{Convergent}, \code{Parallel}}
    \item{\code{Curvature}}{End Curvature: \code{Round}, \code{Medium}, \code{Shallow}}
    \item{\code{Retouched}}{Retouching: \code{Retouched}, \code{Unretouched}}
    \item{\code{Site}}{Site: \code{Castenet A}, \code{Ferrassie H}}
    \item{\code{Freq}}{Number of end scrapers}
  }
}
\details{
The scrapers are grouped on 5 categorical variables into 48 groups. Sackett's analysis employed Chi square and the examination of residuals. Read used the data to illustrate loglinear modelling (1974, 2007). The data come from Tables IV and VIII (pp 373 and 380) in Sackett's original article
}
\source{
Sackett, James R. 1966. Quantitative Analysis of Upper Paleolithic Stone Tools. \emph{American Anthropologist} 68(2): 356--394.
}
\references{
Carlson, David L. 2017. \emph{Quantitative Methods in Archaeology Using R}. Cambridge University Press, pp 72-77, 91-94.

Read, Dwight W. 1974. Some Comments on Typologies in Archaeology and an Outline of a Methodology. \emph{American Antiquity} 39: 216-242.

Read, Dwight W. 2007. \emph{Artifact Classification: A Conceptual and Methodological Approach}. Left Coast Press.
}
\examples{
data(EndScrapers)
xtabs(Freq~Site+Curvature, EndScrapers)
xtabs(Freq~Curvature+Sides+Site, EndScrapers)
}
\keyword{datasets}
