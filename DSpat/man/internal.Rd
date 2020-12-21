\name{Internal}
\alias{Ops.psp}
\alias{rev_val}
\alias{im.clipped}
\alias{AIC.dspat}
\alias{print.dspat}
\alias{vcov.dspat}
\alias{summary.dspat}
\alias{coef.dspat}
\alias{owin.gpc.poly}
\title{Internal DSpat functions}
\description{
Miscellaneous set of functions used in the package.

}
\usage{
\method{AIC}{dspat}(object,...,k)

\method{print}{dspat}(x,...)

\method{summary}{dspat}(object,...)

\method{coef}{dspat}(object,...)

\method{vcov}{dspat}(object,...)

\special{x.psp==y.psp}

\special{x.psp!=y.psp}

rev_val(x,y,val)

im.clipped(x, window)

owin.gpc.poly(window)
}
\arguments{
  \item{x}{for generic functions: a dspat object output from \code{\link{dspat}} with class \code{dspat}, for \code{rev_val}: a vector of x coordinates,
            for \code{im.clipped}: a vector of image values in order defined by \code{spatstat}}
  \item{object}{a dspat object output from \code{\link{dspat}} with class \code{dspat}}
  \item{k}{penalty per parameter in AIC, default is 2}
  \item{x.psp, y.psp}{psp objects}
  \item{y}{vector of y coordinates}
  \item{val}{vector of image values}
  \item{window}{class owin polygonal window}
  \item{...}{additional arguments for generic functions}
  
}
\details{
Internal functions:

 \tabular{ll}{
AIC \tab Computes AIC value\cr
print.dspat \tab print various objects in dspat object\cr
summary.dspat \tab shows summary of ppm model object\cr
coef.dspat \tab provides coefficients of the intensity and detection function\cr
vcov.dspat \tab provides var-cov matrix of coefficients\cr
Ops.psp \tab provides == and != operators for psp objects\cr
\code{rev_val} \tab reverses order of \code{val} such that y increases within increasing x as needed in \code{im}\cr
im.clipped \tab creates image and fills values (\code{val}) into the clipped portion \cr 
           \tab          of the image as defined by window\cr
owin.gpc.poly \tab creates a gpc class poly from first polygon of an owin class\cr
}

}
\author{Jeff Laake}


