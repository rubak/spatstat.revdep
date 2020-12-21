\name{DSpat-package}
\alias{DSpat-package}
\alias{DSpat}
\docType{package}
\title{
Spatial modelling package for distance sampling data
}
\description{
\code{DSpat} uses the tools in \code{spatstat} to provide an analysis of
distance sampling data in a spatial context in which the density surface and
the detection function are estimated simultaneously. The package provides
a fitted density surface and total abundance and measures of precision.  It also
provides simulation capabilities.}
\details{
\tabular{ll}{
Package: \tab DSpat\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2008-04-08\cr
License: \tab GPL version 2 or later\cr
}
Conventional distance sampling (Buckland et al. 2001;2004) uses likelihood theory
for estimation of the detection function based on an assumed uniform distribution
of perpendicular distances within the transects.  An adequate sampling design provides
the basis for the uniform distribution assumption and inference for abundance. No assumption
is made about the spatial distribtion of the object being sampled.

\code{DSpat} provides a full-likelihood framework for simultaneous estimation of the
detection function and abundance based on an inhomogeneous Poisson process.  A
full-likelihood approach has a number of advantages because there is no strict
requirement on the sampling design so it can be used with unequal coverage sampling and
it can provide estimates of the density surface and abundance for any defined sub-area.
Also, by modelling the observed data as a spatial process 'adjustments' to the
strip-width of the transect occur naturally when the transect extends beyond the
area containing objects. Consider sampling a marine environment with a contorted coastline
such as fjords.  In sampling the open ocean, the full transect width can contain objects
but within a fjord the strip is narrowed or clipped in areas where it extends onto land.
This causes difficulties with conventional distance sampling which assumes
a uniform distribution of objects across the entire strip.  Thus, either
a very narrow strip must be used for both areas or the detection function must be
estimated separately for each region and even that can not completely cope with the problem.
This variation in the spatial distribution of objects is handled easily with a spatial
model that simultaneously estimates the detection function and the intensity (density)
of the point process (e.g., animal/plant locations).  The detection function
is estimated as a covariate to explain the intensity of the observed point process as a function of
perpendicular distance from the centerline.  Thus, obviously the potential for
confounding occurs if the pattern of transects is such that pattern of perpendicular distances
is confounded with the spatial pattern of a covariate that determines the true intensity
of the process.  For example, if there was a density gradient with respect to the
coastline and a single-sided transect parallel to the coastline was sampled then perpendicular
distance and distance from the coastline are completely confounded.  However, with a
typical dual-sided transect, the pattern of perpendicular distance is no longer entirely confounded with
the distance from the coastline because perpendicular increases away from the centerline in
both directions.  Thus, confounding would not occur except in the unlikely situation that
intensity (density) varied relative to the coastline in such a fashion that was symmetric
with respect to the centerline of the transect.  With a modicum of care in the design, confounding
between perpendicular distance and spatial covariates can be avoided but the analyst
should always be cognizant of the potential for confounding.

Current Limitations:
  1) assumes no overlap among strips
  2) no handling of cluster size
  3) assumes detection probability on the transect centerline is 1
  4) can only use a detection function of the form \eqn{log(g(x))=h(x)} where \eqn{h(x)} is
     linear in the parameters. For example, \eqn{h(x)=-tau*(distance^2)/2}.  Note that
     any parameter such as \eqn{tau} is not constrained so this does allow for the
     possibility of an increasing detection function.

The first limitation will require some thought and work as we are unaware of any solutions
at present.  If there is overlap, when \code{owin} in \code{spatstat} is called with the \code{poly=transects},
the code will fail.  It is easy to get around this problem to fit the model by
using \code{study.area} as the boundary but the calls to \code{Kinhom} and
\code{\link{lgcp.correction}} will not work properly.  Also, there are some philosophical
and inference issues that need to be considered if sample overlap.  For example, is the
point process fixed during sampling or should the replicate (and overlapping) samples be
considered as independent realizations of the point process.  Even though most designs
do not have overlapping transects in theory, in practice if the line is composed of
contiguous line segments that vary slightly in angle, the transects will overlap when
created from the line segments. Some solution is needed as this is will likely occur
in most real applications.

The latter three limitations can be resolved with the extension of the likelihood and
additional coding in the package.  \code{DSpat} currently uses \code{ppm} in \code{spatstat}
which uses either \code{\link{glm}} or \code{gam} in \code{mgcv} to solve for the MLEs.
We have functions that compute the likelihood and they can be generalized to accomodate these
limitations but they have not been incorporated into the package yet.
  
\code{DSpat} relies heavily on the tools in \code{spatstat} and to a lesser degree the functions in
\code{gpclib}, \code{mgcv} and \code{RandomFields}.  \code{DSpat} provides
aditional functions to cope with analysis and simulation of distance sampling data
(line transect only at this stage). The functions in \code{DSpat} are listed below in
categories with a brief description.

There are a number of concepts that should be understood prior to using this package.
There are 2 coordinate systems that we will use.  The first is the standard x,y Cartesian coordinate
system with x on the horizontal and y on the vertical.  The second which is not used
extensively (yet) is the coordinate system within each line-transect.  A line-transect
is composed of a line (centerline) which has a beginning (x0,y0) and end (x1,y1) and a
rectangular strip with a defined \code{width} which extends \code{width/2} to the left and right
of the centerline.  We use the term line to represent the line and transect for the rectangular
strip (line-transect). The transect has a left-half and right-half defined by the direction of
travel from beginning to the end of the line.  Imagine the line-transect rotated such that
it is vertical with the rotated versions of y0,y1 such that y0<y1 (travelling from south to north).
We define a coordinate system u,v within the line-transect.  The origin for u,v (u=0,v=0)
is the rotated location of the beginning of the line (x0,y0) and u
is equivalent to the standard horizontal x-coordinate with a range of \code{(-width/2,width/2)}
and v is equivalent to the vertical y-coordinate with a range of \code{(0,L)} where \code{L}
is the length of the line \code{L=sqrt((x0-x1)^2+(y0-y1)^2)}.  We use the variable \code{distance} for
the perpendicular distance which is the absolute value of u.

So why have 2 coordinate systems?  \code{spatstat} always works with the x,y coordinate system and it
creates grids and the like with a horizontal-vertical orientation to the grid. In fitting distance
sampling data we want to control the grid resolution relative to the u,v coordinates.  In particular,
we need to use a relatively fine grid in the u direction for estimation of the detection function which
can change quickly over a small scale relative to most covariates that would be used for the
intensity function.  To use the \code{spatstat} code for grids and the like, we rotate the line-transects
and observations to vertical from south to north and create the grid and counting weights in what is now
the u,v coordinate system. Thus, for clarity we use a function argument \code{epsvu} in place of
\code{epsyx} to show that the grid resolution is over u,v and not over x,y, unless the line-transects
are all originally oriented vertically. Currently all line-transects must be rectangular we envision
generalizing this and the u,v coordinate system will be used.

Even though all transects must be rectangular, the surveyed portion of the
transect need not be rectangular.  This is relevant when portions of the transect
extend outside the boundary of the study area (defined region being sampled with the transects).
The study area can be defined by any polygon as defined for class \code{owin} in \code{spatstat}.
Note the restriction that the polygon coordinates must be given in a counter-clockwise direction.
A simple example would be a square region such as 

\code{study.area=owin(xrange=c(0,100),yrange=c(0,100))}

or a square with a missing portion 

\code{study.area=owin(poly=data.frame(x=c(0,40,40,100,100,0),
                                      y=c(0,0,50,50,100,100)))}.

You can examine these by simply typing \code{plot(study.area}.  Regardless, of the study area
shape but depending on the orientation of the transects, portions of the transects can extend
outside of the study area. For example, consider the corners of transects at a 45 degree angle
extending across a rectangular study area.  In many practical applications the width of the
transect is so narrow relative to the dimensions of the study area, that these corners are of no
consequence.  However, in some applications with small scales this can be important.  For example, surveys
of narrow inlets (fjords) or rivers or contorted coastlines or surveys of small areas (see \code{\link{weeds}}).

This is handled in \code{DSpat} by clipping the portion of the rectangular transect that extends
outside of the study area.  The transects are clipped after they are rotated and gridded. This is important
because that ensures the grids \code{spatstat} are positioned the same across all transects.

\bold{Analysis}

\emph{Primary Functions:}

\code{\link{dspat}}  - main function for fitting spatial model to distance sampling data

\code{\link{integrate.intensity}} - computes predicted intensity surface, total abundance and precision with optional correction for over-dispersion

\code{\link{transect.intensity}} - computes predicted and observed counts within each transect in specified
perpendicular distance intervals

\emph{Secondary Functions:}

\code{\link{create.covariate.images}} - create a list of covariate images from a dataframe of covariates.
The list of covariate images is used by \code{\link{LTDataFrame}}.

\code{\link{lines_to_strips}} - from a dataframe of lines this function creates a \code{psp}
object and a list of transect polygons that assumes that lines are the centerlines of strips that have width
as defined in the lines dataframe.

\code{\link{lgcp.correction}} - computes Monte Carlo correction for over-dispersion

\code{\link{LTDataFrame}} - assign covariates to the data (observations) and dummy points

\code{\link{quadscheme.lt}} - constructs a quadrature scheme for \code{ppm} that is
more efficent for line transect samples which are small slices of the study area.  The default
quadrature scheme in \code{spatstat} creates dummy points across the entire study area which
is terribly inefficient.  This function rotates each line to vertical, creates a quadrature
scheme within the line and then "rotates" back to original position to get the proper covariates.
These line-by-line quadratures are then merged into a single quadrature.


\bold{Data preparation and utility}

\code{\link{offset.points}} - this utility function is useful for
most applied data sets in which the position of the observation is specified by the coordinates on the
line that are perpendicular to the object. For a line and its observations, this function converts the object positions on the line and
the perpendicular distance (negative is left) to the coordinates for the location of the
object.  It could be generalized to work with a radial distance and angle which would often be
collected in shipboard work. It is also used from \code{\link{lines_to_strips}} to compute the
vertices of the transect from the lines with a given width.

\code{\link{create.points.by.offset}} - this is a wrapper function  that calls \code{\link{offset.points}} for
each line in a lines dataframe and the corresponding observations in an observations dataframe, and
returns a new observations dataframe with x,y being the true object coordinates.

\code{\link{dist2line}} - this function is the inverse of \code{\link{offset.points}}. It takes the
true coordinates of points and a line and computes the perpendicular distance on the line and
the coordinates on the line.

\code{\link{project2line}} - likewise this is a wrapper function for \code{\link{dist2line}}
that is essentially the inverse of \code{\link{create.points.by.offset}}.


\bold{Internal}

     \code{\link{AIC.dspat}} - computes AIC for the model; only correct if a HPP or IPP process
     
     \code{\link{coef.dspat}} -extracts the coefficients into a list with a vector for intensity coefficients and another
     for detection coefficients.

     \code{\link{print.dspat}} - provides a listing of elements in the dspat object.

     \code{\link{summary.dspat}} - extracts the \code{ppm} object and calls the \code{spatstat} summary function
     for this object.

     \code{\link{vcov.dspat}} - extracts the variance-covariance matrix from the \code{ppm} object.
     
     \code{\link{Ops.psp}} - allows syntax x==y or x!=y where x and y are psp objects.

     \code{\link{rev_val}} - reorders vector for use in \code{im}.

     \code{\link{im.clipped}} - fills in clipped image with vector of values defined
                                over the clipped region.

     \code{\link{owin.gpc.poly}} - converts first polygon in owin class to a gpc polygon.

\bold{Simulation}

\code{\link{create.lines}} - create a systematic grid of parallel lines (with a random start)
across a study area at a specified angle.

\code{\link{sample.points}} - extract observed points from a point process that fall within
the defined set of strips and are randomly detected with a defined detection function.

\code{\link{simCovariates}} - a non-general function for simulating covariates
 in a 100x100 rectangle with discrete habitats and a linear vertical habitat feature. See \code{\link{DSpat.covariates}}.

\code{\link{simDSpat}} - a wrapper function to simulate distance sampling from a rectangular
study area with a specified set of covariates on a grid.  It calls \code{\link{create.lines}},
\code{\link{lines_to_strips}}, \code{\link{simPts}} and \code{\link{sample.points}} and
returns a dataframe of \code{lines} and \code{observations} that can be used with the \code{covariates}
datframe in \code{\link{dspat}} for an analysis.

\code{\link{simPts}} - creates a simulated point process in a study area by calling
\code{RFsimulate} from the package \code{RandomFields} and \code{rpoispp} from
\code{spatstat}.  The intensity process is defined by a \code{covariates} dataframe and a formula
and parameters for the intensity as a function of the covariates.
     

\bold{Example datasets}

An example dataset from the fairy tale simulated world of \code{\link{simCovariates}} can be found in
\code{\link{DSpat.obs}}, \code{\link{DSpat.lines}}, \code{\link{DSpat.covariates}}.  To run an
example analysis with these data, type \code{example(dspat)} or \code{example(DSpat)} to run the same code below.

An example real-world dataset of a devil's claw weed in a farm paddock can be found in
\code{\link{weeds}}, \code{\link{weeds.all}}, \code{\link{weeds.obs}}, \code{\link{weeds.lines}}, \code{\link{weeds.covariates}}.
To run a set of analyses, type \code{example(weeds)}.

}
\author{
 Devin S. Johnson, Jeffrey L. Laake, and Jay M. Ver Hoef

Maintainer: <Jeff.Laake@Noaa.Gov>
}
\references{
Johnson,D.S., Laake, J.L., and Ver Hoef, J.M. (in prep). A model based approach for making ecological inference from distance sampling data.

Buckland, S.T., D.R.Anderson, K.P. Burnham, J.L. Laake, D.L. Borchers, and L. Thomas. 2001.
Introduction to Distance Sampling: Estimating Abundance of Biological Populations. Oxford University Press.

Buckland, S.T., D.R.Anderson, K.P. Burnham, J.L. Laake, D.L. Borchers, and L. Thomas. 2004.
Advanced Distance Sampling. Oxford University Press.

}
\seealso{ \code{\link[spatstat:spatstat-package]{spatstat}} }
\examples{
\donttest{
# get example data
data(DSpat.lines)
data(DSpat.obs)
data(DSpat.covariates)
# Fit model with covariates used to create the data
sim.dspat=dspat(~ river + factor(habitat),
                study.area=owin(xrange=c(0,100), yrange=c(0,100)),
                obs=DSpat.obs,lines=DSpat.lines,covariates=DSpat.covariates,
                epsvu=c(1,.01),width=0.4)
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
}
}

