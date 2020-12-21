\name{RBPottery}
\alias{RBPottery}
\docType{data}
\title{
Romano-British Pottery
}
\description{
Results of chemical analyses of 48 specimens of Romano-British pottery from 5 sites in 3 regions.
}
\usage{data("RBPottery")}
\format{
  A data frame with 48 observations on the following 12 variables.
  \describe{
    \item{\code{ID}}{Sample ID}
    \item{\code{Kiln}}{Kiln: \code{Gloucester}, \code{Llanedeyrn}, \code{Caldicot}, \code{Islands Thorns}, and \code{Ashley Rails}}
    \item{\code{Region}}{Region: \code{Gloucester}, \code{Wales}, and \code{New Forest}}
    \item{\code{Al2O3}}{Percentage aluminum trioxide}
    \item{\code{Fe2O3}}{Percentage Iron trioxide}
    \item{\code{MgO}}{Percentage magnesium oxide}
    \item{\code{CaO}}{Percentage calcium oxide}
    \item{\code{Na2O}}{Percentage sodium oxide}
    \item{\code{K2O}}{Percentage potassium oxide}
    \item{\code{TiO2}}{Percentage titanium dioxide}
    \item{\code{MnO}}{Percentage manganese oxide}
    \item{\code{BaO}}{Percentage barium oxide}
  }
}
\details{
Results of chemical analyses of 48 specimens of Romano-British pottery published by Tubb, et al. (1980). The numbers are the percentage metal oxide. "Kiln" indicates at which kiln site the pottery was found. The kiln sites come from three regions (1=Gloucester, (2=Llanedeyrn, 3=Caldicot), (4=Islands Thorns, 5=Ashley Rails)). The data were scanned from Table 2.2 in Baxter (2003, p. 21) and preserve three probable typographical errors in the original publication. Those errors are the values for TiO2 in line 4 (sample GA4), for MnO in line 35 (sample C13), and for K2O in line 36 (sample C14). Versions of these data are also available as \code{Pottery} in package \code{car}, \code{pottery} in package \code{HSAUR}, and \code{Pottery2} in package \code{heplots}.
}
\source{
Baxter, M. J. 2003. \emph{Statistics in Archaeology}. Arnold.

Tubb, A., A. J. Parker, and G. Nickless. 1980. The Analysis of Romano-British Pottery by Atomic Absorption Spectrophotometry. \emph{Archaeometry} 22: 153-71.
}

\references{
Carlson, David L. 2017. \emph{Quantitative Methods in Archaeology Using R}. Cambridge University Press, pp 247-255, 335-342.
}
\examples{
data(RBPottery)
print(aggregate(RBPottery[, -c(1:3)], list(Region=RBPottery$Region), mean), digits=2)
plot(Na2O~CaO, RBPottery, pch=as.numeric(Region)-1)
legend("topright", levels(RBPottery$Region), title="Region", pch=0:2)
}
\keyword{datasets}
