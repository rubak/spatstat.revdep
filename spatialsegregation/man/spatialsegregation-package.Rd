\name{spatialsegregation-package}
\alias{spatialsegregation-package}
\alias{spatialsegregation}
\docType{package}

\title{Spatial Segregation Measures}

\description{
Collection of measures or summaries of spatial multitype exposure: segregation vs. mingling of different
types of points in a plane.
}


\details{

This is a collection of summaries for multitype spatial point patterns (see package \link{spatstat} for more). 

The package is developed for an article Rajala\&Illian 2010, and provides summaries for detecting simple inter-type 
effects in the pattern.

See the help of the functions for further information.  

Package provides an example dataset object called \code{exposurepps}, documented separately.

Also, the Dixon bivariate test based on contingency tables is available.

\preformatted{
Functions
------------------
segregationFun   -   General calculation function, please use one of the following wrappers:
minglingF        -   Mingling index
shannonF         -   Spatial Shannon index
simpsonF         -   Spatial Simpson index
isarF            -   ISAR function

mingling.index   -   Shortcut for a single value
shannon.index    -   -''-
simpson.index    -   -''-
isar.index       -   -''-

dixon            -   Dixon's 2-type contingency table tests
}

 
}

\author{
Tuomas Rajala
tuomas.rajala@iki.fi
}
\examples{
data(exposurepps)
help(exposurepps)
}
\references{
Graz: The behaviour of the species mingling index \code{$m_{sp}$} in relation to species dominance and dispersion. Eur. J. forest research. 123:87-92, 2004.

Lewandowski, Pommerening: Zur Beschreibung der Waldstruktur - Erwartete und beobachtete Arten-Durchmischung. Forstwiss Centralbl, 116:129-139, 1997.

Rajala, Illian: A family of spatial biodiversity measures based on graphs, Env. Ecol. Stat. 2012

Reardon, O'sullivan: Measures of spatial segregation. Sociological methodology, 34:121-162, 2004.

Shimatani, Kubota: Quantitative assesment of multispecies spatial pattern with high species diversity. Ecological Research, 19, 2004.

Wiegand, Gunatilleke, Gunatilleke, Huth: How individual species structure diversity in tropical forests. PNAS, nov 16, 2007. 
}

\keyword{ package }
