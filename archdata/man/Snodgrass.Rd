\name{Snodgrass}
\alias{Snodgrass}
\docType{data}
\title{
House pits at the Mississippian Snodgrass site in Butler County, Missouri, U.S.A.
}
\description{
Information on the size, location and contents of 91 house pits at the Snodgrass site which was occupied between about CE 1325-1420.
}
\usage{data(Snodgrass)}
\format{
  A data frame with 91 observations on the following 15 variables.
  \describe{
    \item{\code{East}}{East grid location of house in feet (excavation grid system)}
    \item{\code{South}}{East grid location of house in feet (excavation grid system)}
    \item{\code{Length}}{House length in feet}
    \item{\code{Width}}{House width in feet}
    \item{\code{Segment}}{Three areas within the site  \code{1}, \code{2}, \code{3}}
    \item{\code{Inside}}{Location within or outside the "white wall" \code{Inside}, \code{Outside}}
    \item{\code{Area}}{Area in square feet}
    \item{\code{Points}}{Number of projectile points}
    \item{\code{Abraders}}{Number of abraders}
    \item{\code{Discs}}{Number of discs}
    \item{\code{Earplugs}}{Number of earplugs}
    \item{\code{Effigies}}{Number of effigies}
    \item{\code{Ceramics}}{Number of ceramics}
    \item{\code{Total}}{Total Number of artifacts listed above}
    \item{\code{Types}}{Number of kinds of artifacts listed above}
  }
}
\details{
The data from 91 house pits at the Snodgrass site were reported by Price and Giffin in 1979. The layout of the houses follows a grid pattern with the long axis oriented northeast surrounded by a fortification trench. There is also evidence of an interior wall that may have separated the houses inside that wall from those outside the wall. Price and Griffin use differences in house size and artifact composition to suggest that those differences may have reflected rank differences between the occupants of the two areas. That conclusion has been questioned on a number of grounds by Cogswell, et al (2001), but the data are still useful for illustrating a number of quantitative methods. The data come from the appendices except for the house locations which were estimated from the base map in Figure 10 (Price and Griffin 1979).
}
\source{
Price, J. E. and J. B. Griffin. 1979. The Snodgrass Site of the Powers Phase of Southeast Missouri. \emph{Anthropological Papers. Museum of Anthropology, University of Michigan, No.} 66.
}
\references{
Carlson, David L. 2017. \emph{Quantitative Methods in Archaeology Using R}. Cambridge University Press, pp 171-183, 232-242.

Cogswell, J. W., M. J. O'Brien, and D. S. Glover. 2001. The Artifactual Content of Selected House Floors at Turner and Snodgrass. In \emph{Mississippian Community Organization: The Powers Phase in Southeastern Missouri}, edited by M. J. O'Brien, pp 181--229. Kluwer Academic/Plenum.
}
\examples{
data(Snodgrass)
plot(-South~East, Snodgrass, main="Snodgrass Site", pch=as.numeric(Inside)+4, asp=1)
legend("topleft", levels(Snodgrass$Inside), pch=5:6)
boxplot(Area~Inside, Snodgrass)
}
\keyword{datasets}
