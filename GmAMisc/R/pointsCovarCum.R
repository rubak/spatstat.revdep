#' R function to plot the cumulative distribution (and acceptance interval) of the values of a
#' spatial covariate measured at the locations of a point pattern
#'
#' The function allows to test if there is a significant dependence of the input point pattern on a
#' underlying spatial numeric covariate (first-order effect).\cr The function takes as input three
#' datasets: a point patter ('SpatialPointsDataFrame' class), a covariate layer (of 'RasterLayer'
#' class), and (optionally) a polygon feature ('SpatialPolygonsDataFrame' class) representing the
#' study area and exactly matching the extent of the covariate layer. If the latter is not provided,
#' it is internally worked out from the covariate raster and may make the whole function take a
#' while to complete.\cr
#'
#' The function plots the cumulative distribution of the values of the covariate at the locations of
#' the input point pattern, and adds an acceptance interval (with significance level equal to 0.05;
#' sensu Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press
#' 2016, 208) that allows to assess the statistical significance of the observed cumulative
#' distribution. The interval is built by calculating the cumulative distribution of B realizations
#' of a Complete Spatial Random process, and keeping the middle 95percent of those B distributions.
#' B is set by default to 200, but can be increased by the user. The number of random points drawn
#' during each of the B simulations is equal to the number of features of the input point
#' pattern.\cr
#'
#' For an example of the cumulative distribution plot plus acceptance interval, \strong{see} for
#' instance Carrero-Pazos, M. (2018). Density, intensity and clustering patterns in the spatial
#' distribution of Galician megaliths (NW Iberian Peninsula). Archaeological and Anthropological
#' Sciences. https://doi.org/10.1007/s12520-018-0662-2, figs. 4 and 5.\cr
#'
#' @param feature Feature (of point type; 'SpatialPointsDataFrame' class) representing the spatial
#'   point pattern of interest.
#' @param cov.var Numeric covariate (of 'RasterLayer' class).
#' @param studyplot Feature (of polygon type; 'SpatialPolygonsDataFrame' class) representing the
#'   study area and exactly matching the extent of the covariate layer. If NULL, it is worked out
#'   from the covariate layer (may make the whole function take a while to complete).
#' @param B Number of randomized iterations to be used to calculate the acceptance interval (200 by
#'   default).
#' @param cov.var.name Name of the input covariate to be used in the cumulative distribution chart
#'   as label for the x axis (NULL by default).
#' @param oneplot Set to TRUE (default), will plot the charts into a single visualization.
#'
#' @return  The function returns a 2 plots, which can be arranged in just one visualization
#' setting the parameter 'oneplot' to TRUE:\cr
#'
#' -a plot of the point pattern against the underlaying covariate;\cr
#'
#' -a plot of the cumulative distribution of the values of the covariate at the locations of the
#' point patter along with the above-mentioned acceptance interval.\cr
#'
#' @keywords pointsCovarCum
#'
#' @export
#'
#' @importFrom raster rasterToPolygons extract
#'
#' @examples
#' #load the point dataset representing the location of Starbucks shops
#' data(Starbucks)
#'
#' #load the polygon dataset representing the study area
#' data(Massachusetts)
#'
#' #load the raster representing the population density, to be used as covariate
#' data(popdensity)
#'
#' results <- pointsCovarCum(feature=Starbucks, cov.var=popdensity, studyplot=Massachusetts,
#' cov.var.name="population density")
#'
#' @seealso \code{\link{pointsCovarModel}}
#'
pointsCovarCum <- function(feature, cov.var, studyplot=NULL, B=200, cov.var.name=NULL, oneplot=TRUE){
  options(scipen=999)

#if studyplot is NULL, studyplot is worked out from the covariate raster, otherwise do nothing
if(is.null(studyplot)==TRUE){
  #create a copy of the covariate raster
  cov.var.copy <- cov.var
  #substitute the values of the copy with a single value, omitting the NA cells
  cov.var.copy[na.omit(cov.var.copy)] <- 1
  #convert the modified copy to a polygon ignoring the NA cells and dissolving cells with the same
  #attribute (that is all, since we changed all the values to 1) into a single polygon
  studyplot <-rasterToPolygons(cov.var.copy, na.rm=TRUE, dissolve=TRUE)
} else {}

region <- studyplot

#for each point in the input feature dataset extract the value of the cavariate
cov.var.extract <- raster::extract(cov.var, feature)

#calculate the ECDF of the covariate values
cov.var.extract.ecdf <- ecdf(cov.var.extract)

#create a matrix to store the covariate value measured at random locations;
#each column correspond to a random set of points
cov.var.rnd.mtrx <- matrix(nrow=length(feature), ncol=B)

#set the progress bar to be used later on within the loop
pb <- txtProgressBar(min = 0, max = B, style = 3)

  for (i in 1:B){
    #draw a random sample of points within the study region
    rnd <- spsample(region, n=length(feature), type='random')
    #extract the covariate values at the random locations and store them in the matrix (column-wise)
    cov.var.rnd.mtrx[,i] <- raster::extract(cov.var, rnd)
    setTxtProgressBar(pb, i)
  }

# Make a list for the ecdfs
rnd.ecdfs <- list()
for(i in 1:ncol(cov.var.rnd.mtrx)){
  rnd.ecdfs[[i]] <- ecdf(cov.var.rnd.mtrx[,i])
}

#set the limit of the x-axis; NA removal is activated since it can happen that some points (either observed or random) can have no value
#for the covariate; for instance, when there has been an error of some sort in generating the covariate raster
xlim = c(min(min(cov.var.rnd.mtrx, na.rm=TRUE), min(cov.var.extract, na.rm=TRUE)), max(cov.var.extract, na.rm=TRUE))

# We will evaluate the ecdfs on a grid of 1000 points between
# the x limits
xs <- seq(xlim[1], xlim[2], length.out = 1000)

# This actually gets those evaluations and puts them into a matrix
out <- lapply(seq_along(rnd.ecdfs), function(i){rnd.ecdfs[[i]](xs)})
tmp <- do.call(rbind, out)

# Get the .025 and .975 quantile for each column
# at this point each column is a fixed 'x' and the rows
# are the different ecdfs applied to that
lower <- apply(tmp, 2, quantile, probs = .025)
upper <- apply(tmp, 2, quantile, probs = .975)

#conditionally set the layout in just one visualization
if(oneplot==TRUE){
  m <- rbind(c(1,2))
  layout(m)
}

#plot the map showing the feature against the covariate
raster::plot(cov.var,
             main="Map of the point pattern against the covariate",
             cex.main=0.95,
             axes=TRUE)

#add the point feature
raster::plot(feature,
             pch=20,
             add=TRUE)

#set the title for the ecdf plot
maintitle <- paste0("Cumulative distribution of the covariate values at the points' location \n(acceptance interval based on ", B, " randomized iterations)")

#set the ecdf plot xlab title
xlab_string <- ifelse(is.null(cov.var.name)==TRUE,"", cov.var.name)

# Plot the original data
# plot the ECDF of the first random dataset
graphics::plot(cov.var.extract.ecdf,
               verticals=TRUE,
               do.points=FALSE,
               col="white",
               xlab=xlab_string,
               ylab="Proportion",
               main=maintitle,
               cex.main=0.90,
               xlim= xlim)

# Add in the quantiles
graphics::polygon(c(xs,rev(xs)), c(upper, rev(lower)), col = "#DBDBDB88", border = NA)

graphics::plot(cov.var.extract.ecdf,
               verticals=TRUE,
               do.points=FALSE,
               add=TRUE)
}
