\name{CCmaps}
\alias{CCmaps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Conditioned choropleth maps}
\description{
Conditioned choropleth maps permit the conditioning of a map of a variable on the values of one or two other variables coded as factors or shingles. This function uses \link[sp]{spplot} after constructing multiple subsets of the variable of interest defined by the intervals given by the conditioning variables.
}

\usage{
CCmaps(obj, zcol = NULL, cvar = NULL, cvar.names = NULL, ..., names.attr,
 scales = list(draw = FALSE), xlab = NULL, ylab = NULL,
 aspect = mapasp(obj, xlim, ylim), sp.layout = NULL, xlim = bbox(obj)[1, ],
 ylim = bbox(obj)[2, ])
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ object of class \link[sp]{SpatialPolygonsDataFrame}}
  \item{zcol}{single variable name as string}
  \item{cvar}{a list of one or two conditioning variables, which should be of class factor or shingle}
  \item{cvar.names}{names for conditioning variables, if not given, the names of the variables in the \code{cvar} list}
  \item{\dots}{other arguments passed to \link[sp]{spplot} and \link[lattice]{levelplot}}
  \item{names.attr}{names to use in panel, if different from zcol names}
  \item{scales}{scales argument to be passed to Lattice plots; use
\code{list(draw = TRUE)} to draw axes scales}
  \item{xlab}{label for x-axis}
  \item{ylab}{label for y-axis}
  \item{aspect}{aspect ratio for spatial axes; defaults to "iso" (one unit
on the x-axis equals one unit on the y-axis) but may be set to more
suitable values if the data are e.g. if coordinates are latitude/longitude}
  \item{sp.layout}{NULL or list; see \link[sp]{spplot}}
  \item{xlim}{numeric; x-axis limits}
  \item{ylim}{numeric; y-axis limits}
}

\value{
 The function returns a \link[sp]{SpatialPolygonsDataFrame} object with the \code{zcol} variable and the partitions of the \code{cvars} list variables invisibly.
}
\references{
Carr D, Wallin J, Carr D (2000) Two new templates for epidemiology applications: linked micromap plots and conditioned choropleth maps. \emph{Statistics in Medicine} 19(17-18): 2521-2538
Carr D, White D, MacEachren A (2005) Conditioned choropleth maps and hypothesis generation. \emph{Annals of the Association of American Geographers} 95(1): 32-53
Friendly M (2007) A.-M. Guerry's Moral Statistics of France: challenges for multivariable spatial analysis. \emph{Statistical Science} 22(3): 368-399
}
\author{ Roger Bivand }

\seealso{\code{\link[sp]{spplot}}}

\examples{
nc.sids <- readShapeSpatial(system.file("shapes/sids.shp",
 package="maptools")[1], IDvar="FIPSNO",
 proj4string=CRS("+proj=longlat +ellps=clrk66"))
nc.sids$ft.SID74 <- sqrt(1000)*(sqrt(nc.sids$SID74/nc.sids$BIR74) +
 sqrt((nc.sids$SID74+1)/nc.sids$BIR74))
nc.sids$ft.NWBIR74 <- sqrt(1000)*(sqrt(nc.sids$NWBIR74/nc.sids$BIR74) +
 sqrt((nc.sids$NWBIR74+1)/nc.sids$BIR74))
library(lattice)
sh_nw4 <- equal.count(nc.sids$ft.NWBIR74, number=4, overlap=1/5)
CCmaps(nc.sids, "ft.SID74", list("Nonwhite_births"=sh_nw4),
 col.regions=colorRampPalette(c("yellow1", "brown3"))(20),
 main="Transformed SIDS rates 1974-8")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

