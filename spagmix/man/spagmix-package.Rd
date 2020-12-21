\name{spagmix-package}
\Rdversion{1.1}
\alias{spagmix-package}
\alias{spagmix}
\docType{package}
\title{
The spagmix Package: Artificial Spatial and Spatiotemporal Densities on Bounded Windows
}
\description{
Provides functions to design synthetic spatial and spatiotemporal densities and relative risk functions based mainly on Gaussian mixture distributions, and simulate iid data therefrom.
}
\details{
\tabular{ll}{
Package: \tab spagmix\cr
Version: \tab 0.3-2\cr
Date: \tab 2020-06-24\cr
License: \tab GPL (>= 2) \cr
}
Appraisal of existing, refined, and new statistical methods for the analysis of spatial and spatiotemporal point pattern data usually involves numeric experimentation. Motivated by relevant problems in nonparametric density estimation (see e.g. Wand & Jones, 1995), \code{spagmix} (``spatial Gaussian mixtures'') provides some simple utilities for designing heterogeneous density and density-ratio or \emph{relative risk} (Bithell 1990, 1991; Kelsall & Diggle, 1995) functions in space and space-time (see Fernando & Hazelton, 2014 for the latter). The package is also capable of producing realisations of (possibly inhomogeneous) spatial log-Gaussian Cox process intensities (\enc{Møller}{Moller} et al., 1998; see also Davies & Hazelton, 2013).

Additionally, the package contains functions for simulating datasets given these scenarios. For examples of how these kinds of synthetic functions have been used in simulation studies in various publications, see for example Clark & Lawson, 2004; Davies & Hazelton, 2010; Davies, 2013a,b; Davies & Hazelton, 2013; Fernando et al., 2014; Davies et al., 2016; Davies et al., 2018a; and Davies & Lawson, 2019.

We have designed the objects of \code{spagmix} to use and be compatible with standard object classes of the \code{\link[spatstat:spatstat-package]{spatstat}} (Baddeley & Turner, 2005; Baddeley et al., 2015) and \code{\link[sparr:sparr-package]{sparr}} (Davies et al., 2018b) packages. The content of \code{spagmix} can be broken up as follows:\cr

\emph{Artificial polygonal windows}\cr
Some pre-made synthetic spatial windows; these are all single closed polygons as objects of class \code{\link[spatstat.geom]{owin}} and are lazy-loaded:\cr
\code{\link{bx}}, \code{\link{heart}}, \code{\link{shp1}}, \code{\link{shp2}}, \code{\link{star}}, \code{\link{toywin}}\cr


\emph{Spatial scenarios}\cr
\code{\link{sgmix}} is used to create spatial (2D) Gaussian mixture distributions on a bounded subset of the plane.\cr
\code{\link{rgmix}} also creates 2D Gaussian mixture densities, but does so by stochastic generation of the contributing bumps.\cr
\code{\link{rrmix}} creates Gaussian mixture relative risk scenarios based on a supplied control density (see e.g. Davies & Hazelton, 2010).
\code{\link{lgcpmix}} generates a spatial log-Gaussian Cox process intensity in space, given a deterministic intensity function and residual correlation governed by a stochastic realisation of a Gaussian field with a specified covariance structure.

\emph{Spatiotemporal scenarios}\cr
\code{\link{stgmix}} is used to create spatiotemporal (3D) Gaussian mixture densities on a bounded subset of the plane and a single closed interval in time.\cr
\code{\link{stkey}} is used to create spatiotemporal densities by pixel-wise interpolation of multiple spatial image `keyframes'.
\code{\link{rrstmix}} is a spatiotemporal version of \code{\link{rrmix}}, used to create artificial spatiotemporal relative risk functions. Note the control density may be purely spatial, representing a distribution `at-risk' points that does not change over time (Fernando & Hazelton, 2014).\cr


\emph{Data generation}\cr
To generate purely spatial data for a single spatial density, the user is directed to \code{\link[spatstat.core]{rpoint}} of the \code{\link[spatstat:spatstat-package]{spatstat}} package or \code{\link[sparr]{rimpoly}} of the \code{\link[sparr:sparr-package]{sparr}} package.\cr
\code{\link{rpoispoly}} is a wrapper of \code{\link[sparr]{rimpoly}}, and is used to generate realisations of Poisson point processes in space, given an intensity function.
\code{\link{rrpoint}} is a wrapper of \code{\link[sparr]{rimpoly}}, and is used to generate iid datasets based on a synthetic spatial relative risk surface object.\cr
\code{\link{rstpoint}} is a 3D rejection algorithm for sampling iid data from a supplied spatiotemporal density.\cr
\code{\link{rrstpoint}} is a wrapper of \code{\link{rstpoint}} to generate iid datasets from a synthetic spatiotemporal relative risk surface object.\cr

\emph{Miscellaneous}\cr
\code{\link{plot.stim}} is an \code{S3} plotting method for spatiotemporal density objects.\cr
\code{\link{stintegral}} computes the 3D integral of a spatiotemporal density object.\cr
\code{\link{unify.owin}} is a wrapper for \code{\link[spatstat.geom]{affine}} that transforms any spatial \code{\link[spatstat.geom]{owin}} to fall inside the unit square.
}

\section{Dependencies/Imports}{
The \code{spagmix} package depends upon \code{\link[spatstat:spatstat-package]{spatstat}} (Baddeley & Turner, 2005; Baddeley et al., 2015) and imports from \code{RandomFields} (Schlather et al., 2015), \code{abind} (Plate & Heiberger, 2016), \code{sparr} (Davies et al., 2018b), and mvtnorm (Genz et al., 2018). We also highly recommend the \code{rgl} package (Adler et al., 2018) which can be used to create interactive plots of spatiotemporal data.
}

\author{
A.K. Redmond and T.M. Davies\cr
\emph{Dept. of Mathematics & Statistics, University of Otago, Dunedin, New Zealand}\cr

Maintainer: T.M.D. \email{tdavies@maths.otago.ac.nz}\cr
}


\references{
Adler, D., Murdoch, D. and others (2018), rgl: 3D Visualization Using OpenGL, R package version 0.99.16 https://CRAN.R-project.org/package=rgl\cr\cr
Baddeley, A., Rubak, E. and Turner, R. (2015), \emph{Spatial Point Patterns: Methodology and Applications with R}, Chapman and Hall/CRC Press, London.\cr\cr
Baddeley, A. and Turner, R. (2005), Spatstat: an R package for analyzing spatial point patterns, \emph{Journal of Statistical Software}, \bold{12}(6), 1-42.\cr\cr
Bithell, J.F. (1990), An application of density estimation to geographical epidemiology, \emph{Statistics in Medicine}, \bold{9}, 691-701.\cr\cr
Bithell, J.F. (1991), Estimation of relative risk function,. \emph{Statistics in Medicine}, \bold{10}, 1745-1751.\cr\cr
Clark, A.B. and Lawson, A.B. (2004), An evaluation of non-parametric relative risk estimators for disease maps, \emph{Computational Statistics & Data Analysis}, \bold{47}, 63-78.\cr\cr
Davies, T.M. (2013a), Jointly optimal bandwidth selection for the planar kernel-smoothed density-ratio, \emph{Spatial and Spatio-temporal Epidemiology}, \bold{5}, 51-65.\cr\cr
Davies, T.M. (2013b), Scaling oversmoothing factors for kernel estimation of spatial relative risk, \emph{Epidemiological Methods}, \bold{2}(1), 67-83.\cr\cr
Davies, T.M. and Hazelton, M.L. (2010), Adaptive kernel estimation of spatial relative risk, \emph{Statistics in Medicine}, \bold{29}(23), 2423-2437.\cr\cr
Davies, T.M. and Hazelton, M.L. (2013), Assessing minimum contrast parameter estimation for spatial and spatiotemporal log-Gaussian Cox processes, \emph{Statistica Neerlandica}, \bold{67}(4), 355-389.\cr\cr
Davies, T.M., Jones, K. and Hazelton, M.L. (2016), Symmetric adaptive smoothing regimens for estimation of the spatial relative risk function, \emph{Computational Statistics & Data Analysis}, \bold{101}, 12-28.\cr\cr
Davies, T.M. and Lawson, A.B. (2019), An evaluation of likelihood-based bandwidth selectors for spatial and spatiotemporal kernel estimates, \emph{Journal of Statistical Computation and Simulation}, \bold{89} 1131-1152.\cr\cr
Davies, T.M., Flynn, C.R. and Hazelton, M.L. (2018a), On the utility of asymptotic bandwidth selectors for spatially adaptive kernel density estimation, \emph{Statistics & Probability Letters}, \bold{138}, 75-81.\cr\cr
Davies, T.M., Marshall, J.C. and Hazelton, M.L. (2018b), Tutorial on kernel estimation of continuous spatial and spatiotemporal relative risk, \emph{Statistics in Medicine}, \bold{37}(7), 1191-1221.\cr\cr
Fernando, W.T.P.S., Ganesalingam, S. and Hazelton, M.L. (2014), A comparison of estimators of the geographical relative risk function, \emph{Journal of Statistical Computation and Simulation}, \bold{84}(7), 1471-1485.\cr\cr
Fernando, W.T.P.S. and Hazelton, M.L. (2014), Generalizing the spatial relative risk function, \emph{Spatial and Spatio-temporal Epidemiology}, \bold{8}, 1-10.\cr\cr
Genz, A., Bretz, F., Miwa, T., Mi, X., Leisch, F., Scheipl, F. and Hothorn, T. (2018), mvtnorm: Multivariate Normal and t Distributions, R package version 1.0-8. URL http://CRAN.R-project.org/package=mvtnorm \cr\cr
Kelsall, J.E. and Diggle, P.J. (1995), Kernel estimation of relative risk, \emph{Bernoulli}, \bold{1}, 3-16.\cr\cr
\enc{Møller}{Moller}, J., Syversveen, A.R. and Waagepetersen, R.P. (1998), Log-Gaussian Cox processes, \emph{Scandinavian Journal of Statistics}, \bold{25}(3) 451--482.\cr\cr
Plate, T. and Heiberger, R. (2016), abind: Combine Multidimensional Arrays, R package version 1.4-5. https://CRAN.R-project.org/package=abind \cr\cr
Schlather, M., Malinowski, A., Menck, P.J., Oesting, M. and Strokorb, K. (2015) Analysis, Simulation and Prediction of Multivariate Random Fields with Package RandomFields, \emph{Journal of Statistical Software}, \bold{63}(8) 1--25.\cr\cr
Wand, M.P. and Jones, C.M. (1995), \emph{Kernel Smoothing}, Chapman & Hall, London.
}

\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
%\examples{

%}
