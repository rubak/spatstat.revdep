\name{rrstmix}
\alias{rrstmix}
\alias{rrstim}
\title{
Spatiotemporal relative risk surface generation
}
\description{
Generates an appropriately scaled spatiotemporal (trivariate) relative risk surface using a supplied control density and \eqn{N} Gaussian-style hotspots.
}
\usage{
rrstmix(g, rhotspots, rsds, rweights, rbase = 1, log = TRUE,
        tlim = NULL, tres = NULL)
}
\arguments{
  \item{g}{
The control density as a \code{\link{stim}}, \code{\link[sparr:spattemp.density]{stden}}, or \code{\link[spatstat]{im}} object; this will be internally rescaled to integrate to 1 if it does not already do so. When a \code{\link{stim}} or \code{\link[sparr:spattemp.density]{stden}} object, the resolution and domain of the final result will be the same as this. When this argument is passed an object of class \code{\link[spatstat]{im}}, the function assumes a static (unchanging) control density over time (see Fernando and Hazelton, 2014), and the user must additionally specify \code{tlim} and \code{tres}.
}
  \item{rhotspots}{
A \eqn{3 \times N} matrix specifying the spatiotemporal coordinates of the \eqn{N} peaks and troughs in the relative risk density. The three entries down each column will be respectively interpreted as \eqn{x}-coord., \eqn{y}-coord., and time-coordinate of each Gaussian bump.
}
  \item{rsds}{
A \eqn{3 \times N} strictly positive numeric matrix specifying the standard deviations along each axis of each of the \eqn{N} bumps, the ordering of the components in each column is the same as \code{rhotspots}.
}
  \item{rweights}{
  A vector of length \eqn{N} giving relative weightings for each peak (positive weight) or trough (negative).
}
\item{rbase}{
The base level of the relative risk surface (default is 1). The peaks and troughs will be added or subtracted from this base level prior to normalisation.
}
  \item{log}{
  A logical value. If \code{TRUE} (default), the relative risk surface is returned logged.
}
  \item{tlim}{
Only used if \code{g} is a pixel \code{\link[spatstat]{im}}age object. A vector of length 2 giving the boundaries of the time interval on which the relative risk surface will be defined.
}
  \item{tres}{
Only used if \code{g} is a pixel \code{\link[spatstat]{im}}age object. The resolution along the temporal axis of the final result.
}
}

\details{
This function is the spatiotemporal (trivariate) equivalent of \code{\link{rrmix}}. See `Details' in the documentation for that function for more information.
}

\references{
Fernando, W.T.P.S. and Hazelton, M.L. (2014), Generalizing the spatial relative risk function, \emph{Spatial and Spatio-temporal Epidemiology}, \bold{8}, 1-10.
}

\value{
An oject of class \code{rrstim}. This is a list with the following components:
\item{f}{An object of class \code{\link{stim}} giving the case density.}
\item{g}{A copy of the object passed to the argument \code{g}, possibly renormalised to integrate to 1 if this was necessary. If \code{g} was originally an \code{\link[spatstat]{im}}, this will be converted to an object of class \code{\link{stim}}.}
\item{r}{An object of class \code{\link{stim}} giving the (log) relative risk surface.}
}
\author{A.K. Redmond and T.M. Davies}



\examples{
\donttest{
# time-varying control density
gg1 <- stgmix(mean=matrix(c(2,1,3,0,-1,5),nrow=3),
             vcv=array(c(1,0,0,0,1,0,0,0,1,2,0,0,0,1,0,0,0,2),dim=c(3,3,2)),
             window=shp2,tlim=c(0,6))
rsk1 <- rrstmix(g=gg1,rhotspots=matrix(c(-2,0,2,1,2,5.5),nrow=3),
               rsds=sqrt(cbind(rep(1.5,3),rep(0.25,3))),rweights=c(-0.5,5))
plot(rsk1$g,sleep=0.1,fix.range=TRUE)
plot(rsk1$f,sleep=0.1,fix.range=TRUE)
plot(rsk1$r,sleep=0.1,fix.range=TRUE)

# time-constant control density
set.seed(321)
gg2 <- rgmix(7,window=shp2)
rsk2 <- rrstmix(g=gg2,rhotspots=matrix(c(-1,-1,2,2.5,0,5),nrow=3),
                rsds=sqrt(cbind(rep(0.75,3),c(0.05,0.01,0.5))),
                rweights=c(-0.4,7),tlim=c(0,6),tres=64)
plot(rsk2$g,sleep=0.1,fix.range=TRUE)
plot(rsk2$f,sleep=0.1,fix.range=TRUE)
plot(rsk2$r,sleep=0.1,fix.range=TRUE)
}
}

