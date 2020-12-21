\name{OxfordPots}
\alias{OxfordPots}
\docType{data}
\title{
Distribution of Late Romano-British Oxford Pottery
}
\description{
Percentages of Late Romano-British Oxford Pottery on 30 sites
}
\usage{data("OxfordPots")}
\format{
  A data frame with 30 observations on the following 7 variables.
  \describe{
    \item{\code{Place}}{Site name}
    \item{\code{OxfordPct}}{Percentage of Oxford pottery}
    \item{\code{OxfordDst}}{Distance to Oxford in miles}
    \item{\code{NewForestPct}}{Percentage of New Forest pottery}
    \item{\code{NewForestDst}}{Distance to New Forest}
    \item{\code{WalledArea}}{Acreage of walled towns}
    \item{\code{WaterTrans}}{Availability of a water transportation link, 1=probable presence}
  }
}
\details{
In several publications Ian Hodder analyzed the spatial distribution of Late Romano-British pottery produced at Oxford as evidence of trade and marketing patterns. These data come from the article by Fulford and Hodder (1974). In addition to the percentage of Oxford pottery and the distance to Oxford for 30 sites, data on New Forest pottery was included and information on walled town size and the availability of water transportation.
}
\source{
Fulford, M. and I. Hodder. 1974. A Regression Analysis of Some Late Romano-British Pottery: A Case Study. \emph{Oxoniensia} 39: 26-33.
}
\references{
Hodder, I. 1974. A Regression Analysis of Some Trade and Marketing Patterns. \emph{World Archaeology} 6: 172-189.

Hodder, I. and C. Orton. 1976. \emph{Spatial Analysis in Archaeology}, pp 117-119.
}
\examples{
data(OxfordPots)
# Construct Fulford and Hodder's Figure 3
Oxford.lm1 <- lm(log(OxfordPct)~OxfordDst, OxfordPots, subset=WaterTrans==0)
Oxford.lm2 <- lm(log(OxfordPct)~OxfordDst, OxfordPots, subset=WaterTrans==1)
plot(log(OxfordPct)~OxfordDst, OxfordPots, xlim=c(0, 160), yaxt="n", ylim=c(0, 3.25),
     ylab="Percentage of Oxford Pottery", xlab="Distance (miles)", 
     pch=c(1, 16)[WaterTrans+1], cex=1.5, lwd=2)
# Add log y-axis
axis(2, log(c(1, 5, 10, 20)), c(1, 5, 10, 20), las=1)
abline(Oxford.lm1, lwd=2)
abline(Oxford.lm2, lwd=2)
}
\keyword{datasets}
