\name{Handaxes}
\alias{Handaxes}
\docType{data}
\title{
Lower Paleolithic handaxes from Furze Platt, Maidenhead, Berkshire, England
}
\description{
Handaxes from the Furze Platt site stored at the Royal Ontario Museum.
}
\usage{data(Handaxes)}
\format{
  A data frame with 600 observations on the following 8 variables.
  \describe{
    \item{\code{Catalog}}{Specimen catalog number}
    \item{\code{L}}{Maximum Length}
    \item{\code{L1}}{Distance from the butt to the location of the maximum breadth measured along the length dimension}
    \item{\code{B}}{Maximum breadth}
    \item{\code{B1}}{Breadth measured at 1/5 of the length from the tip. Measured perpendicular to the length}
    \item{\code{B2}}{Breadth measured at 1/5 of the length from the butt. Measured perpendicular to the length}
    \item{\code{T}}{Maximum thickness, not necessarily measured at the maximum breadth}
    \item{\code{T1}}{Thickness measured at B1}
  }
}
\details{
The data consist of measurements on 600 handaxes from the Furze Platt site stored at the Royal Ontario Museum that were measured by William Fox. The measurements follow the system used by Derek Roe (Roe 1964, 1968, 1981). Fox's measurements were digitized by Tony Baker and uploaded to his website.
}
\source{
Fox, William and Tony Baker. 2006. Dimensions of 600 Acheulean Handaxes from Furze Platt, Maidenhead, Berkshire, England. Online at \url{http://ele.net/acheulean/FPatROM.htm}.
}
\references{
Baker, Tony. 2006. The Acheulean Handaxe. Online article available at \url{http://ele.net/acheulean/handaxe.htm}. Accessed Oct 14, 2013.

Carlson, David L. 2017. \emph{Quantitative Methods in Archaeology Using R}. Cambridge University Press, pp 221-231, 269-277, 321-330.

Fox, William. 1969. An Analysis of the R. O. M. Collection, Lower Paleolithic Implements: Furze Platt, Maidenhead, Berkshire, England. Unpublished paper and notes in the possession of the author.

Roe, D. A. 1964. The British Lower and Middle Paleolithic: some problems, methods of study and preliminary results. \emph{Proceedings of the Prehistoric Society} 30: 245--267.

Roe, D. A. 1968. British Lower and Middle Paleolithic Handaxe groups. \emph{Proceedings of the Prehistoric Society} 34: 1--82.

Roe, D. A. 1981. \emph{The Lower and Middle Paleolithic Periods in Britain}. Routledge.
}
\examples{
data(Handaxes)
summary(Handaxes)
}
\keyword{datasets}
