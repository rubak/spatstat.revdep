\name{rpoispoly}
\alias{rpoispoly}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate a Poisson point pattern in a polygonal window
}
\description{
Generates a single realisation of a spatial Poisson point process based on a pixel \code{\link[spatstat]{im}}age and a polygonal \code{\link[spatstat]{owin}}.
}
\usage{
rpoispoly(z, w = NULL, correction = 1.1, maxpass = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{A pixel image of class \code{\link[spatstat]{im}} defining the spatial intensity function of the points. The number of points generated, \eqn{n}, will be found as a randomly generated Poisson variate with mean parameter equal to the integral of \code{z}.}
\item{w}{A polygonal window of class \code{\link[spatstat]{owin}}. See `Details'.}

\item{correction}{An adjustment to the number of points generated at the initial pass of the internal loop in an effort to minimise the total number of passes required to reach \eqn{n} points. See `Details' and `Warning'.}

\item{maxpass}{The maximum number of passes allowed before the function exits. If this is reached before \eqn{n} points are found that fall within \code{w}, a warning is issued.}
}
\details{
This is a wrapper function for \code{\link[sparr]{rimpoly}} that operates much like \code{\link[spatstat]{rpoispp}}, but with artificial corrections at the edges of boundary pixels. This allows the user to generate a realisation of a 2D Poisson point process using a supplied pixel \code{\link[spatstat]{im}}age as the spatial intensity function, but return the result with a polygonal \code{\link[spatstat]{owin}} instead of a binary image mask.

Let \eqn{n} be a randomly generated integer from a Poisson distribution with mean given by the integral of the intensity function \code{z}. When the user specifies their own polygonal window in \code{w}, a \code{while} loop is called and repeated as many
times as necessary (up to \code{maxpass} times) to find \eqn{n} points inside \code{w} (when \code{w = NULL}, then the union of the pixels of \code{z}
is used, obtained via \code{as.polygonal(Window(z))}). The loop is necessary because the standard behaviour of \code{\link[spatstat]{rpoispp}} can (and often does)
yield points that sit in corners of pixels which lie outside a corresponding irregular polygon \code{w}.

The \code{correction} argument is used to determine how many points are generated initially,
which will be \code{ceiling(correction*n)}; to minimise the number of required passes over the loop this is by default set to give a number slightly higher than the requested \eqn{n}.

An error is thrown if \code{Window(z)} and \code{w} do not overlap.
}
\value{
An object of class \code{\link[spatstat]{ppp}} containing the Poisson-generated points, defined with the polygonal \code{\link[spatstat]{owin}}, \code{w}.
}

\references{
Diggle, P.J. (2014) \emph{Statistical Analysis of Spatial and Spatiotemporal Point Patterns}, 3rd Ed, Chapman & Hall, Boca Raton, USA.
}

\author{
T.M. Davies
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\section{Warning}{
Note that this is an artificial correction that forces the Poisson-generated number of \eqn{n} points to be found inside \emph{any} supplied polygon \code{w} (even if \code{w} only partially covers the domain of \code{z}). As such, this function only makes sense in terms of the theory of a Poisson point process if the polygon \code{w} corresponds exactly to the pixellised intensity. For practical intents and purposes, it therefore must be assumed in using this function that a supplied polygon \code{w} is/was the original basis for the discretisation into the pixel image for the purposes of producing the intensity \code{z}, and hence that any adverse effects arising from imposing \code{w} as the window of the final result are negligible. See `Examples'.
}

\seealso{
\code{\link[spatstat]{rpoint}}, \code{\link[sparr]{rimpoly}}, \code{\link[spatstat]{rpoispp}}
}

\examples{
mn <- cbind(c(0.25,0.8),c(0.31,0.82),c(0.43,0.64),c(0.63,0.62),c(0.21,0.26))
v1 <- matrix(c(0.0023,-0.0009,-0.0009,0.002),2)
v2 <- matrix(c(0.0016,0.0015,0.0015,0.004),2)
v3 <- matrix(c(0.0007,0.0004,0.0004,0.0011),2)
v4 <- matrix(c(0.0023,-0.0041,-0.0041,0.0099),2)
v5 <- matrix(c(0.0013,0.0011,0.0011,0.0014),2)
vr <- array(NA,dim=c(2,2,5))
for(i in 1:5) vr[,,i] <- get(paste("v",i,sep=""))
intens <- sgmix(mean=mn,vcv=vr,window=toywin,p0=0.1,int=500)

aa <- rpoispp(intens) # Default spatstat function
bb <- rpoispoly(intens) # No polygon supplied; just uses pixel union
cc <- rpoispoly(intens,w=toywin) # Original irregular polygon

par(mfrow=c(2,2))
plot(intens,log=TRUE)
plot(aa,main=paste("aa\nn =",npoints(aa)))
plot(bb,main=paste("bb\nn =",npoints(bb)))
plot(cc,main=paste("cc\nn =",npoints(cc)))
}
