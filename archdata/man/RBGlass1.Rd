\encoding{UTF-8}
\name{RBGlass1}
\alias{RBGlass1}
\docType{data}
\title{
Romano-British Glass, Major and Minor Elements
}
\description{
The concentrations for 11 major and minor elements in 105 Romano-British waste glass specimens from two furnace sites (Leicester and Mancetter).
}
\usage{data("RBGlass1")}
\format{
  A data frame with 105 observations on the following 12 variables.
  \describe{
    \item{\code{Site}}{a factor with levels \code{Leicester} and \code{Mancetter}}
    \item{\code{Al}}{Percentage Aluminum}
    \item{\code{Fe}}{Percentage Iron}
    \item{\code{Mg}}{Percentage Magnesium}
    \item{\code{Ca}}{Percentage Calcium}
    \item{\code{Na}}{Percentage Sodium}
    \item{\code{K}}{Percentage Potassium}
    \item{\code{Ti}}{Percentage Titanium}
    \item{\code{P}}{Percentage Phosphorus}
    \item{\code{Mn}}{Percentage Manganese}
    \item{\code{Sb}}{Percentage Antinmony}
    \item{\code{Pb}}{Percentage Lead}
  }
}
\details{
The concentrations for 11 major and minor elements in 105 Romano-British waste glass specimens from two furnace sites (Leicester and Mancetter) come from Caroline Jackson's Ph. D. thesis at Bradford University. The data here were scanned from from Baxter (1994) Table A1. Measurements are percentage for each element.
}
\source{
Baxter, M. J. 1994. \emph{Exploratory Multivariate Analysis in Archaeology}. Edinburgh University Press.

Jackson, C. M. 1992. A Compositional Analysis of Roman and Early Post-Roman Glass and Glass Working Waste from Selected British Sites Towards an Understanding of the Technology of Glass-Making Through Analysis by Inductively-Coupled Plasma Spectrometry. Unpublished PhD thesis. Bradford University (BL: D214554).
}
\references{
Baxter, M. J., Cool H.E.M., Heyworth M.P. and Jackson, C.M. 1995. Compositional Variability in Colourless Roman Vessel Glass. \emph{Archaeometry} 37(1), 129-141.

Baxter, M. J., Cool, H. E. M. and Jackson, C. M. (2005). Further Studies in the Compositional Variability of Colourless Romano-British Glass.  \emph{Archaeometry}  47, 47-68.

Carlson, David L. 2017. \emph{Quantitative Methods in Archaeology Using R}. Cambridge University Press, pp 245-247, 256-261.

Jackson, C M, J R Hunter, S E Warren, and H E M Cool. 1991. The Analysis of Blue-Green Glass and Glassy Waste from Two Romano-British Glass Working Sites. In \emph{Archaeometry 1990}, edited by E. Pernicka and G. A. Wagner, pp 295-304. Birkhäuser Verlag.
}
\examples{
data(RBGlass1)
RBGlass1.pca <- prcomp(RBGlass1[, -1], scale.=TRUE)
biplot(RBGlass1.pca, xlabs=abbreviate(RBGlass1$Site, 1), cex=.75)
}
\keyword{datasets}
