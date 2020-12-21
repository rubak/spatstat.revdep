\name{rrpoint}
\alias{rrpoint}
\alias{rrstpoint}
\title{
Generate random case/control points in space or space-time
}
\description{
Generates a pair of random, independent point patterns corresponding to a case density and a control density, for relative risk analyses.
}
\usage{
rrpoint(n, r, W = NULL, correction = 1.1, maxpass = 50)
rrstpoint(n, r, W = NULL, correction = 1.5, maxpass = 50)
}
\arguments{
  \item{n}{
The number of points to be generated. This must be a numeric vector of length 2 giving the number of points to generate for the case and control densities respectively. Alternatively a single number can be supplied; then the same number of points is generated for both densities.
}
  \item{r}{
The relative risk surface object containing the definitions of the case and control probability densities: an object of class \code{\link{rrim}} or \code{\link[sparr:risk]{rrs}} for \code{rrpoint}, or an object of class \code{\link{rrstim}} or \code{\link[sparr:spattemp.risk]{rrst}} for \code{rrstpoint}.
}
  \item{W}{
The polygonal \code{\link[spatstat.geom]{owin}} defining the spatial window on which the density is defined. If \code{NULL}, this will be set to the \code{\link[spatstat.geom]{as.polygonal}} version of the pixel images stored in \code{r}. See `Details'.
}
  \item{correction}{
An adjustment to the number of points generated at the initial pass of the internal loop in an effort to minimise the total number of passes required to reach \code{n} points.
}
  \item{maxpass}{
The maximum number of passes allowed before the function exits. If this is reached before \code{n} points are found with respect to the spatial or spatiotemporal domains of \code{r}, a warning is issued.
}
}
\details{
These functions randomly generate a pair of independent spatial or spatiotemporal point patterns of \code{n} points based on the case and control density functions stored in \code{r}. At any given pass for each density, \code{n} * \code{correction} points are generated and rejection sampling is used to accept some of the points; this is repeated until the required number of points is found.

The argument \code{W} is optional, but is useful when the user wants the spatial window of the resulting point pattern to be a corresponding irregular polygon, as opposed to being based on the boundary of a binary image mask (which, when the pixel \code{\link[spatstat.geom]{im}}ages in \code{r} are converted to a polygon directly, gives jagged edges based on the union of the pixels).
}
\value{
A list with two components, \code{cases} and \code{controls}, each of which is an object of class \code{\link[spatstat.geom]{ppp}} containing the \code{n} generated points. for spatiotemporal densities, the \code{\link[spatstat.geom]{marks}} of the object will contain the correspondingly generated observation times.
}
\author{
T.M. Davies
}
\examples{
# Using 'rrim' object:
set.seed(1)
gg <- rgmix(3,window=toywin,S=matrix(c(0.08^2,0,0,0.1^2),nrow=2),p0=0.2)
rho <- rrmix(g=gg,
             rhotspots=cbind(c(0.8,0.3),c(0.4,0.4),c(0.6,0.5),c(0.3,0.5)),
             rsds=c(0.005,0.025,0.01,0.025),
             rweights=c(3,2,10,5)*10)

rho.sample <- rrpoint(n=c(400,800),r=rho,W=toywin)

par(mfrow=c(2,2))
plot(rho$g,main="control density")
plot(rho$f,main="case density")
plot(rho$r,main="log relative risk surface")
plot(rho.sample$controls,main="sample data")
points(rho.sample$cases,col=2)
legend("topright",col=2:1,legend=c("cases","controls"),pch=1)


\donttest{
# Using 'rrs' object:
require("sparr")
data(pbc)
pbccas <- split(pbc)$case
pbccon <- split(pbc)$control
h0 <- OS(pbc,nstar="geometric")
f <- bivariate.density(pbccas,h0=h0,hp=2,adapt=TRUE,pilot.density=pbccas,
                       edge="diggle",davies.baddeley=0.05,verbose=FALSE)
g <- bivariate.density(pbccon,h0=h0,hp=2,adapt=TRUE,pilot.density=pbccas,
                       edge="diggle",davies.baddeley=0.05,verbose=FALSE)
pbcrr <- risk(f,g,tolerate=TRUE,verbose=FALSE)

pbcrr.pt <- rrpoint(n=1000,r=pbcrr)

par(mfrow=c(1,3))
plot(pbcrr)
plot(pbcrr.pt$cases)
plot(pbcrr.pt$controls)


# Using 'rrstim' object:
set.seed(321)
gg <- rgmix(7,window=shp2)
rsk <- rrstmix(g=gg,rhotspots=matrix(c(-1,-1,2,2.5,0,5),nrow=3),
               rsds=sqrt(cbind(rep(0.75,3),c(0.05,0.01,0.5))),
               rweights=c(-0.4,7),tlim=c(0,6),tres=64)
plot(rsk$r,fix.range=TRUE)

rsk.pt <- rrstpoint(1000,r=rsk,W=shp2)

par(mfrow=c(1,2))
plot(rsk.pt$cases)
plot(rsk.pt$controls)


# Using 'rrst' object:
require("sparr")
data(fmd)
fmdcas <- fmd$cases
fmdcon <- fmd$controls

f <- spattemp.density(fmdcas,h=6,lambda=8)
g <- bivariate.density(fmdcon,h0=6)
rho <- spattemp.risk(f,g)

rho.pt <- rrstpoint(1000,r=rho)

par(mfrow=c(1,2))
plot(rho.pt$cases)
plot(rho.pt$controls)
}
}
