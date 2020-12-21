\name{dspat}
\alias{dspat}
\title{Fits spatial model to distance sampling data}
\description{
Creates a dspat object by fitting model represented by formula to observations
along line transects in a study area with covariates defined for a grid over the
entire study area.
}
\usage{
dspat(int.formula=~1, det.formula=~1, study.area, obs, lines, covariates,
       epsvu=c(1,.01), width=NULL, use.gam=FALSE, show.warnings=FALSE,
       nclass=NULL)
}
\arguments{
  \item{int.formula}{formula for model of the point process intensity}
  \item{det.formula}{formula for interaction with distance in the detection process}
  \item{study.area}{owin class for study area}
  \item{obs}{dataframe of observations}
  \item{lines}{dataframe of lines}
  \item{covariates}{dataframe of covariates on a grid in the study area}
  \item{epsvu}{vector of height of pixels(y) and width of pixels(x)}
  \item{width}{full transect width; only needed if it is not specified in lines.df}
  \item{use.gam}{if TRUE uses gam instead of glm for fitting; if formula contains s() use.gam will be set TRUE by default}
  \item{show.warnings}{if TRUE, show the warnings created in building the quadrature.}
  \item{nclass}{number of distance classes for expected/observed counts.}
}
\value{  list of class "dspat" with elements
  \item{model}{output object from ppm}
  \item{lines.psp}{psp line segment process for center lines}
  \item{transects}{list of dataframes specifying polygonal transects}
  \item{covariate.im}{list of covariate images (class im)}
  \item{study.area}{owin class of study area}
  \item{use.gam}{TRUE if gam used and FALSE otherwise}
}
\details{
\preformatted{
     covariates has following structure
            x   - x coordinate of midpoint of grid cell
            y   - y coordinate of midpoint of grid cell
            ... - any number of covariate

      the data are ordered by column from left to right and
      from bottom to top such that y changes first from smallest
      to largest. Below are matrices showing y,x and their order
               3,1 3,2 3,3   3 6 9
               2,1 2,2 2,3   2 5 8
               1,1 1,2 1,3   1 4 7}
               
The default for the intensity formula (\code{int.formula}) is ~1,
a homogeneous Poisson process. Note that what is actually fitted is ~-1+constant
where constant is 1 everywhere.  This is done to avoid a glitch in \code{vcov.ppm}.
The detection formula (\code{det.formula}) is expressed as a formula that
interacts with \eqn{ I(-distance^2/2)}.  The default of ~1 is a detection function that
is constant everywhere.  If you use ~-1, it will drop distance which assumes a
strip transect with perfect detection within the strip. The variables contained
in \code{int.formula} must be all contained within \code{covariates} because they
need to be defined across the entire study area. The variables contained in
\code{det.formula} can be in \code{covariates} or in \code{lines} because for
prediction of the intensity, distance is set to zero, so these covariates need not
be known across the entire survey area.

The value of \code{epsvu[2]} is adjusted such that it is an even multiple of \code{width/2}
so that the grid points are evenly distributed in the direction of perpendicular distance.
}
\author{Jeff Laake; Devin Johnson }
\seealso{ \code{\link{quadscheme.lt}},\code{\link{LTDataFrame}}}
\examples{

# get example data
data(DSpat.lines)
data(DSpat.obs)
data(DSpat.covariates)
# Fit model with covariates used to create the data
sim.dspat=dspat(~ river + factor(habitat),
                study.area=owin(xrange=c(0,100), yrange=c(0,100)),
                obs=DSpat.obs,lines=DSpat.lines,covariates=DSpat.covariates,
                epsvu=c(4,.1),width=0.4)
\donttest{
# Print
sim.dspat
# Summarize results
summary(sim.dspat)
# Extract coefficients
coef.intensity <- coef(sim.dspat)$intensity
coef.detection <- coef(sim.dspat)$detection
# Extract variance-covariance matrix (inverse information matrix)
J.inv <- vcov(sim.dspat)
# Compute AIC
AIC(sim.dspat)
# Visualize intensity (no. animals per area) and estimate abundance
mu.B <- integrate.intensity(sim.dspat,dimyx=100)
cat('Abundance =       ', round(mu.B$abundance,0), "\n")
dev.new()
plot(mu.B$lambda, col=gray(1-c(1:100)/120), main='Estimated Intensity')
plot(sim.dspat$model$Q$data,add=TRUE)
plot(owin(poly=sim.dspat$transect),add=TRUE)
plot(sim.dspat$lines.psp,lty=2,add=TRUE)

# Compute se and confidence interval for abundance without over-dispersion
mu.B <- integrate.intensity(sim.dspat,se=TRUE,dimyx=100)
cat("Standard Error =  ", round(mu.B$precision$se,0), "\n",
    "95 Percent Conf. Int. =   (", round(mu.B$precision$lcl.95,0), ',',
           round(mu.B$precision$ucl.95,0), ")", "\n")
# Compute se and confidence interval for abundance with over-dispersion estimate
dev.new()
# The rest of the example has been put into a function to speed up package checking; remove # to run
# to run type do.dspat()
do.spat=function()
{
mu.B <- integrate.intensity(sim.dspat,se=TRUE,od=TRUE,reps=30,dimyx=100)
cat("Standard Error (corrected) =     ", round(mu.B$precision.od$se,0), "\n",
    "95 Percent Conf. Int. (corrected) =      (", round(mu.B$precision.od$lcl.95,0),
              ",", round(mu.B$precision.od$ucl.95,0), ")", "\n")
# Fit model with smooth of x and y
sim.dspat=dspat(~ s(x) + s(y),study.area=owin(xrange=c(0,100), yrange=c(0,100)),
                obs=DSpat.obs,lines=DSpat.lines,covariates=DSpat.covariates,
                epsvu=c(1,.01),width=0.4)
AIC(sim.dspat)
# Visualize intensity (no. animals per area) and estimate abundance
mu.B <- integrate.intensity(sim.dspat,dimyx=100,se=TRUE)
cat('Abundance =       ', round(mu.B$abundance,0), "\n")
cat("Standard Error =     ", round(mu.B$precision$se,0), "\n",
    "95 Percent Conf. Int. =      (", round(mu.B$precision$lcl.95,0),
              ",", round(mu.B$precision$ucl.95,0), ")", "\n")
dev.new()
plot(mu.B$lambda, col=gray(1-c(1:100)/120), main='Estimated Intensity')
plot(sim.dspat$model$Q$data,add=TRUE)
plot(owin(poly=sim.dspat$transect),add=TRUE)
plot(sim.dspat$lines.psp,lty=2,add=TRUE)
#
# Fit model with smooth of x and y with interaction
#
sim.dspat=dspat(~ s(x,y),study.area=owin(xrange=c(0,100), yrange=c(0,100)),
                obs=DSpat.obs,lines=DSpat.lines,covariates=DSpat.covariates,
                epsvu=c(1,.01),width=0.4)
AIC(sim.dspat)
# Visualize intensity (no. animals per area) and estimate abundance
mu.B <- integrate.intensity(sim.dspat,dimyx=100,se=TRUE)
cat('Abundance =       ', round(mu.B$abundance,0), "\n")
cat("Standard Error =     ", round(mu.B$precision$se,0), "\n",
    "95 Percent Conf. Int. =      (", round(mu.B$precision$lcl.95,0),
              ",", round(mu.B$precision$ucl.95,0), ")", "\n")
dev.new()
plot(mu.B$lambda, col=gray(1-c(1:100)/120), main='Estimated Intensity')
plot(sim.dspat$model$Q$data,add=TRUE)
plot(owin(poly=sim.dspat$transect),add=TRUE)
plot(sim.dspat$lines.psp,lty=2,add=TRUE)
}
}
}

