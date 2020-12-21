\encoding{latin1}
\name{figuras}
\alias{fig1}
\alias{fig2}
\alias{fig3}
\alias{figuras}
\docType{data}
\title{Artificial point data.}
\description{
  The three different point patterns in the figure 3.1 of De la Cruz (2008) 
}
\usage{
data(fig1)
data(fig2)
data(fig3)
}
\format{
  A data frame with 87 observations on the following 2 variables.
  \describe{
    \item{\code{x}}{x coordinate}
    \item{\code{y}}{y coordinate}
  }
}
\references{
  De la Cruz, M. 2008. \enc{Métodos para analizar datos puntuales. 
  En: \emph{Introducción al Análisis Espacial de Datos en Ecología y Ciencias Ambientales: Métodos y Aplicaciones } 
  (eds. Maestre, F. T., Escudero, A. y Bonet, A.), pp 76-127. Asociación Española de Ecología Terrestre, Universidad Rey Juan Carlos
  y Caja de Ahorros del Mediterráneo, Madrid.}{Metodos para analizar datos puntuales. 
  En: \emph{Introduccion al Analisis Espacial de Datos en Ecologia y Ciencias Ambientales: Metodos y Aplicaciones } 
  (eds. Maestre, F. T., Escudero, A. y Bonet, A.), pp 76-127. Asociacion Espanola de Ecologia Terrestre, Universidad Rey Juan Carlos
  y Caja de Ahorros del Mediterraneo, Madrid.}
}
\examples{
\dontrun{
data(fig1)

data(fig2)

data(fig3)

# transform to ppp format of spatstat with function haz.ppp:

fig1.ppp <- haz.ppp(fig1)

fig2.ppp <- haz.ppp(fig2)

fig3.ppp <- haz.ppp(fig3)

#Analyses as in Fig.3.2 of De la Cruz (2008). First, compute function K:

cosa1 <- Kest(fig1.ppp)

# Plot different estimators. 
# Fig. 3.2a: 

par("mar"=par("mar")+c(0,1,0,0))

plot(cosa1, col=c(1,0,0,1), lwd=c(2,2,2,2), lty=c(1,1,1,2), 
	 main="")

# Fig. 3.2b:

plot(cosa1, sqrt(./pi)-r~r, col=c(1,0,0,1), lwd=c(2,2,2,2), 
         lty=c(1,1,1,2), main="", ylab="L(r)")

# Fig. 3.2c:

plot(cosa1, .-(pi*r^2)~r, col=c(1,0,0,1), lwd=c(2,2,2,2), 
         lty=c(1,1,1,2), main="", ylab=expression(K(r)-pi*r^2))

# Fig. 3.2d:

plot(cosa1,(./(pi*r^2))-1~r, col=c(1,0,0,1), lwd=c(2,2,2,2), 
	 lty=c(1,1,1,2), main="",
         ylab=expression((K(r)/pi*r^2)-1))

## Analyses as in fig. 3.7 of De la Cruz (2008).
## First, compute function K and pointwise envelopes:

cosa1.env <- envelope(fig1.ppp, Kest)

cosa2.env <- envelope(fig2.ppp, Kest)

cosa3.env <- envelope(fig3.ppp, Kest)

## Plot function L with pointwise envelopes:

plot(cosa1.env,sqrt(./pi)-r~r, lwd=c(1,1,2,2), 
	 lty=c(1,1,3,3), col=c(1,1,1,1), xlab="r", 
         ylab="L(r)", main="", ylim=c(-2,2))

## Add simultaneous envelopes of Ripley (+-1.68 *sqrt(A)/N):

abline(h=1.68*sqrt(fig1.ppp$w$area)/fig1.ppp$n, 
         lty=2, lwd=2)

abline(h=-1.68*sqrt(fig1.ppp$w$area)/fig1.ppp$n, 
         lty=2, lwd=2)

## Plot function L with pointwise envelopes:

plot(cosa2.env,sqrt(./pi)-r~r, lwd=c(1,1,2,2), 
         lty=c(1,1,3,3), col=c(1,1,1,1), xlab="r", 
         ylab="L(r)", main="")

## Add simultaneous envelopes of Ripley:

abline(h=1.68*sqrt(fig2.ppp$w$area)/fig2.ppp$n, 
         lty=2, lwd=2)

abline(h=-1.68*sqrt(fig2.ppp$w$area)/fig2.ppp$n, 
         lty=2, lwd=2)

## Plot function L with pointwise envelopes:

plot(cosa3.env,sqrt(./pi)-r~r, lwd=c(1,1,2,2), 
         lty=c(1,1,3,3), col=c(1,1,1,1), xlab="r", 
         ylab="L(r)", main="")

## Add simultaneous envelopes of Ripley:

abline(h=1.68*sqrt(fig3.ppp$w$area)/fig3.ppp$n, 
         lty=2, lwd=2)

abline(h=-1.68*sqrt(fig3.ppp$w$area)/fig3.ppp$n, 
         lty=2, lwd=2)
}
}
\keyword{datasets}
