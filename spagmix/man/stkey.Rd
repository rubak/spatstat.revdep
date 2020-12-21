\name{stkey}
\alias{stkey}
\title{
Spatiotemporal density generation via keyframe interpolation
}
\description{
Uses the supplied spatial pixel images and scalings to linearly interpolate the behaviour of the function over time, creating a trivariate density function in space-time.
}
\usage{
stkey(start, stop, tlim, kf = NULL, tres = 64,
      kftimes = NULL, fscales = NULL, window = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{start}{
The spatial pixel \code{\link[spatstat]{im}}age corresponding to the spatial density at start of the time interval. May be unnormalised, the function internally rescales all supplied spatial images to integrate to 1.
}
  \item{stop}{
The pixel \code{\link[spatstat]{im}}age for the end of the time interval. Must be compatible with \code{start}, in that it is defined over the same spatial domain and is of identical resolution.
}
  \item{tlim}{
A numeric vector of length 2 representing the temporal window i.e. the time interval over which the interpolation takes place.
}
  \item{kf}{
A \code{\link[spatstat]{solist}} of the pixel \code{\link[spatstat]{im}}ages of the keyframes between \code{start} and \code{stop}. If supplied, each image must be compatible with \code{start} and \code{stop}. If unsupplied, the resulting interpolation is performed only on \code{start} and \code{stop}
}
  \item{tres}{
The resolution of the resulting array in the time dimension (default is 64).
}
  \item{kftimes}{
A vector of times that position the interim keyframes in \code{kf} between \code{tlim[1]} and \code{tlim[2]}. Ignored if \code{kf = NULL}. If unsupplied (\code{NULL}), but \code{kf} is present, the function simply positions the images of \code{kf} at evenly spaced time points in the \code{tlim} interval.
}
  \item{fscales}{
A numeric vector of unnormalised, relative point-intensity scales. This may be provided either as of \code{length(kf)} + 2, so the point intensities assinged to each frame in the order \code{c(start,}<entries of \code{kf}>\code{,stop)}, or of length \code{tres}. If unsupplied, each spatial frame is simply given equal weight. See `Details'.
}
  \item{window}{
An object of class \code{\link[spatstat]{owin}} giving the polygonal spatial observational window on which the density is defined. If \code{NULL}, the polygon is simply obtained from the union of pixel values in the supplied images.
}
}
\details{
This function interpolates in a pixel-wise fashion between the \code{\link[spatstat]{im}} objects supplied as \code{start} and \code{stop} (and \code{kf} if supplied), placing them as keyframes at the times in \code{tlim} (for \code{start} and \code{stop}) and \code{kftimes} (for the members of \code{kf}). The final result is rescaled such that its total integrated volume over the defined spatiotemporal domain is 1, yeilding a trivariate density function.

If \code{fscale} is a vector of length \code{tres}, each element will correspond to the \bold{relative} overall scaling of one of the resulting interpolated pixel images. If it is of length \code{length(kf) + 2}, the scales will correspond to \code{start}, each keyframe in \code{kf} and \code{stop} in that order. The values in this argument are only interpretable in a relative sense: for example, with a single keyframe suppled to \code{kf} (in addition to the required \code{start} and \code{stop}), then \code{fscales = c(0.5,1,0.5)} has exactly the same effect on the final result as \code{fscales = c(1,2,1)}, and is interpreted as yielding a point density that reaches twice the concentration at the time of the supplied keyframe relative to the \code{start} and \code{stop} margins. Supplying \code{fscale} as a vector of length \code{tres} thus allows finer control over the relative point density over time, such as for the incorporation of harmonic seasonal variation.
}


\value{
Like \code{\link{stgmix}}, an object of class \code{\link{stim}} giving the trivariate density. This is a list with six components:
\item{a}{The \code{xr} \eqn{x} \code{yr} \eqn{x} \code{tres} array of the specified density (where \code{xr} and \code{yr} are used here to denote the spatial resolution along the x- and y-axes; this is governed in \code{stkey} by the images initally supplied to \code{start} and \code{stop}).}
\item{v}{A pixel \code{\link[spatstat]{im}}age version of \code{a}, provided as a \code{\link[spatstat]{solist}} of length \code{tres}, with each member being the spatial \code{\link[spatstat]{im}}age slice of the 3D density at each of the time-coordinate values.}
\item{xcol}{Grid coordinates in the spatial x-axis (corresponds to each spatial \code{\link[spatstat]{im}}age in \code{v}).}
\item{yrow}{Grid coordinates in the spatial y-axis (corresponds to each spatial \code{\link[spatstat]{im}}age in \code{v}).}
\item{tlay}{Grid coordinates in the temporal axis (corresponds to the order of the spatial \code{\link[spatstat]{im}}ages in \code{v}).}
\item{W}{A copy of \code{window}, the spatial \code{\link[spatstat]{owin}} upon which the density is defined.}
}

\author{
A.K. Redmond and T.M. Davies
}

\examples{
\donttest{
mn <- matrix(c(0,0,1,2,0.5,-1),nrow=2)
vr <- array(c(0.2,0,0,2,1,0,0,1,1,0.3,0.3,0.5),dim=c(2,2,3))
im1 <- sgmix(mn,c(1,2,1),shp1,p=c(0.4,0.3,0.3))
im2 <- sgmix(matrix(c(-3,0,0,-2,-1,2),nrow=2),c(3,1,1),shp1,p=c(0.4,0.3,0.3))
im3 <- sgmix(mn,vr,shp1,p0=0.1)

kf1 <- stkey(start=im1,stop=im2,tlim=c(5,20),window=shp1)
plot(kf1)

kf2 <- stkey(start=im1,stop=im1,tlim=c(0,15),kf=solist(im1,im1),kftimes=c(2,8),
             fscale=c(1,2,1.5,1),window=shp1)
plot(kf2,fix.range=TRUE)

kf3 <- stkey(start=im1,stop=im2,tlim=c(0,20),kf=solist(im1,im2),kftimes=c(8,12),
             fscale=c(1,2,2,1),window=shp1)
plot(kf3,fix.range=TRUE)

ff <- c(sin((1:64)/3)+1.5)
plot(ff,type="l")
kf4 <- stkey(start=im1,stop=im2,kf=solist(im3),kftimes=25,tlim=c(0,50),fscale=ff,window=shp1)
plot(kf4,fix.range=TRUE)
}
}
