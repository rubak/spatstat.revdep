\name{lgcpmix}
\alias{lgcpmix}

\title{
Generate a spatial log-Gaussian Cox process intensity
}
\description{
Generate a realisation of a (possibly inhomogeneous) log-Gaussian Cox process (LGCP) spatial intensity function with an identifiable mean structure.
}
\usage{
lgcpmix(lambda, covmodel = "exp", covpars = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambda}{
A pixel \code{\link[spatstat]{im}}age giving the deterministic spatial intensity as the mean structure of the process. The generated Gaussian field will match the dimensions, resolution and domain of this object.
}
  \item{covmodel}{
A character string giving the short name of a spatial covariance model available in the \code{\link[RandomFields]{RandomFields}} package (Schlather et al., 2015). The code searches for a function of this name after adding the prefix \code{"RM"}, piggybacking on the internals of \code{\link[spatstat]{rLGCP}}. See `Details'.
}
  \item{covpars}{
A named list of values for the parameters required by the chosen \code{covmodel}. This will normally include \code{var} and \code{scale}. Both variance and scale parameters default to 1 if this argument is left \code{NULL}.
}
}
\details{
This function allows the user to generate a spatial intensity function \eqn{\Gamma} of the form
\deqn{
\Gamma(x) = \lambda(x)\exp[Y(x)]
}
for \eqn{x \in W}, where \eqn{\lambda(x)} (passed to \code{lambda}) is the deterministic spatial intensity over the spatial domain \eqn{W}, and \eqn{Y(x)} is a Gaussian random field on \eqn{W}. This Gaussian field is defined with a particular spatial covariance function (specified in \code{covmodel}) with variance and scale parameters \eqn{\sigma^2} and \eqn{\phi} respectively, as well as any additionally required parameter values (all specified in \code{covpars}).

The mean parameter \eqn{\mu} of the Gaussian field \eqn{Y} is internally fixed at \eqn{-\sigma^2/2}; negative half the variance. This is for identifiability of the mean structure, forcing \eqn{E[Y(x)] = 1} for all \eqn{x \in W} (see theoretical properties in \enc{Møller}{Moller} et al., 1998). This means the deterministic intensity function \eqn{\lambda(x)} is solely responsible for describing fixed heterogeneity in spatial intensity over \eqn{W}, with the randomly generated Gaussian field left to describe residual stochastic spatial correlation. This presents a highly flexible class of model, even with stationarity and isotropy of the Gaussian field itself, and is intuitively sensible in a variety of applications. See Diggle et al. (2005) and Davies & Hazelton (2013) for example.

As such, the pixel \code{\link[spatstat]{im}}age supplied to \code{lambda} as \eqn{\lambda(x)} must be non-negatively-valued and yield a finite integral. The choice of covariance model and correspondingly required parameters as well as actual simulation of the Gaussian field is deferred to functionality in the \code{\link[RandomFields]{RandomFields}} package; see \code{\link[RandomFields]{RMmodel}} for possible choices. For example, requesting \code{covmodel = "exp"} (default) will search for the \code{\link[RandomFields]{RandomFields}} function \code{\link[RandomFields]{RMexp}} and imposes an exponential covariance structure on the generated field whereby \eqn{Cov(u) = \sigma^2\exp(-u/\phi)} for the Euclidean distance between any two spatial locations \eqn{u}.

To generate a subsequent dataset, use e.g. \code{\link[spatstat]{rpoispp}} or \code{\link{rpoispoly}}.
}

\value{
A pixel \code{\link[spatstat]{im}}age giving the generated intensity function, comprised of the product of \code{lambda} (fixed, and unchanging in repeated calls to this function) and the exponentiated Gaussian field (with expected value 1, this is stochastic and varies in repeated calls).
}

\references{
Davies, T.M. and Hazelton, M.L. (2013), Assessing minimum contrast parameter estimation for spatial and spatiotemporal log-Gaussian Cox processes, \emph{Statistica Neerlandica}, \bold{67}(4) 355--389.

Diggle, P.J., Rowlingson, B. and Su, T. (2005), Point process methodology for on-line spatio-temporal disease surveillance, \emph{Environmetrics}, \bold{16} 423--434.

\enc{Møller}{Moller}, J., Syversveen, A.R. and Waagepetersen, R.P. (1998), Log-Gaussian Cox processes, \emph{Scandinavian Journal of Statistics}, \bold{25}(3) 451--482.

Schlather, M., Malinowski, A., Menck, P.J., Oesting, M. and Strokorb, K. (2015) Analysis, Simulation and Prediction of Multivariate Random Fields with Package RandomFields, \emph{Journal of Statistical Software}, \bold{63}(8) 1--25.
}

\author{
T.M. Davies, based partially on code written for the \code{\link[spatstat]{rLGCP}} function by A. Jalilian, R. Waagepetersen, A. Baddeley, R. Turner and E. Rubak.
}


\seealso{
\code{\link[spatstat]{rLGCP}}, \code{\link[spatstat]{rpoispp}}, \code{\link{rpoispoly}}
}

\examples{
## Homogeneous example ##

# Create constant intensity image integrating to 500

homog <- as.im(as.mask(toywin))
homog <- homog/integral(homog)*500


# Corresponding LGCP realisations using exponential covariance structure
par(mfrow=c(2,2),mar=rep(1.5,4))
for(i in 1:4){
  temp <- lgcpmix(homog,covmod="exp",covpars=list(var=1,scale=0.2))
  plot(temp,main=paste("Realisation",i),log=TRUE)
}



## Inhomogeneous examples ##

# Create deterministic trend

mn <- cbind(c(0.25,0.8),c(0.31,0.82),c(0.43,0.64),c(0.63,0.62),c(0.21,0.26))
v1 <- matrix(c(0.0023,-0.0009,-0.0009,0.002),2)
v2 <- matrix(c(0.0016,0.0015,0.0015,0.004),2)
v3 <- matrix(c(0.0007,0.0004,0.0004,0.0011),2)
v4 <- matrix(c(0.0023,-0.0041,-0.0041,0.0099),2)
v5 <- matrix(c(0.0013,0.0011,0.0011,0.0014),2)
vr <- array(NA,dim=c(2,2,5))
for(i in 1:5) vr[,,i] <- get(paste("v",i,sep=""))
intens <- sgmix(mean=mn,vcv=vr,window=toywin,p0=0.1,int=500)


# Two realisations (identical calls to function), exponential covariance structure

r1exp <- lgcpmix(lambda=intens,covmodel="exp",covpars=list(var=2,scale=0.05))
r2exp <- lgcpmix(lambda=intens,covmodel="exp",covpars=list(var=2,scale=0.05))


# Two more realisations, Matern covariance with smoothness 1

r1mat <- lgcpmix(lambda=intens,covmodel="matern",covpars=list(var=2,scale=0.05,nu=1))
r2mat <- lgcpmix(lambda=intens,covmodel="matern",covpars=list(var=2,scale=0.05,nu=1))


# Plot everything, including 'intens' alone (no correlation)

par(mar=rep(2,4))
layout(matrix(c(1,2,4,1,3,5),3))
plot(intens,main="intens alone",log=TRUE)
plot(r1exp,main="realisation 1\nexponential covar",log=TRUE)
plot(r2exp,main="realisation 2\nexponential covar",log=TRUE)
plot(r1mat,main="realisation 1\nMatern covar",log=TRUE)
plot(r2mat,main="realisation 2\nMatern covar",log=TRUE)


# Plot example datasets
dint <- rpoispoly(intens,w=toywin)
d1exp <- rpoispoly(r1exp,w=toywin)
d2exp <- rpoispoly(r2exp,w=toywin)
d1mat <- rpoispoly(r1mat,w=toywin)
d2mat <- rpoispoly(r2mat,w=toywin)

par(mar=rep(2,4))
layout(matrix(c(1,2,4,1,3,5),3))
plot(dint,main="intens alone",log=TRUE)
plot(d1exp,main="realisation 1\nexponential covar",log=TRUE)
plot(d2exp,main="realisation 2\nexponential covar",log=TRUE)
plot(d1mat,main="realisation 1\nMatern covar",log=TRUE)
plot(d2mat,main="realisation 2\nMatern covar",log=TRUE)
}

