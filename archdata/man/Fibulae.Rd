\encoding{UTF-8}
\name{Fibulae}
\alias{Fibulae}
\docType{data}
\title{
Bronze La Tene fibulae from Munsingen, Switzerland
}
\description{
The La Tène fibulae from the Iron Age cemetery of Münsingen near Berne, Switzerland (100 - 500 BCE) described by F. R. Hodson (1968).
}
\usage{data("Fibulae")}
\format{
  A data frame with 30 observations on the following 16 variables.
  \describe{
    \item{\code{Grave}}{Grave number}
    \item{\code{Mno}}{Museum number}
    \item{\code{FL}}{Foot Length}
    \item{\code{BH}}{Bow Height}
    \item{\code{BFA}}{Bow Front Angle}
    \item{\code{FA}}{Foot Angle}
    \item{\code{CD}}{Coil Diameter}
    \item{\code{BRA}}{Bow Rear Angle}
    \item{\code{ED}}{Element Diameter}
    \item{\code{FEL}}{Foot Extension Length}
    \item{\code{C}}{Catchplate}
    \item{\code{BW}}{Bow Width}
    \item{\code{BT}}{Bow Thickness}
    \item{\code{FEW}}{Foot Extension Width}
    \item{\code{Coils}}{Number of Coils}
    \item{\code{Length}}{Total Length}
  }
}
\details{
The La Tène fibulae from the Iron Age cemetery of Münsingen near Berne, Switzerland were reported by F. R. Hodson (1968). They were featured in several papers by Hodson over the years and used to illustrate a variety of multivariate statistical techniques. The data here were taken from Doran and Hodson (1975), Table 9.1. These are the raw  measurements including 5 missing values in foot extension thickness and 1 in foot extension length.
}
\source{
Doran, J. E. and F. R. Hodson. 1975. \emph{Mathematics and Computers in Archaeology}. Harvard University Press, Cambridge, Massachusetts.

Hodson, F. R. 1968. \emph{The La Tène Cemetery at Műnsingen-Rain}. Stampfli, Berne.
}
\references{
Carlson, David L. 2017. \emph{Quantitative Methods in Archaeology Using R}. Cambridge University Press, pp 88-91, 95-99, 103-109, 127-129, 132-138, 162-169.

Hodson, F. R., P. H. A. Sneath, J. E. Doran. 1966. Some Experiments in the Numerical Analysis of Archaeological Data.  \emph{Biometrika} 53: 311-324.

Hodson, F. R. 1969. Searching for Structure within Multivariate Archaeological Data. \emph{World Archaeology} 1: 90-105.

Hodson, F. R. 1970. Cluster Analysis and Archaeology: some New Developments and Applications. \emph{World Archaeology} 1: 299-320.

Hodson, F. R. 1971. Numerical Typology and Prehistoric Archaeology. In \emph{Mathematics int eh Archaeological and Historical Sciences}, edited by F. R. Hodson, D. G. Kendall and P. Tautu, pp 30-45. Edinburgh University Press, Edinburgh.

Sneath, P. H. A. 1968. Goodness of Intuitive Arrangements into Time Trends Based on Complex Pattern. \emph{Systematic Zoology} 17: 256-260.

}
\examples{
data(Fibulae)
t(sapply(Fibulae[, 3:16], quantile, na.rm=TRUE))
plot(density(Fibulae$Length, bw="SJ"), main="Kernel Density Plot of Length")
}
\keyword{datasets}
