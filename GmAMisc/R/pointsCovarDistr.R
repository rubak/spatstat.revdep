#' R function to plot the frequency distribution of the average value of a
#' spatial covariate measured at randomized locations
#'
#' The function allows to test if there is a significant dependence of the input point pattern on a
#' underlying numeric covariate (first-order effect).\cr The function takes as input three
#' datasets: a point patter ('SpatialPointsDataFrame' class), a covariate layer (of 'RasterLayer'
#' class), and (optionally) a polygon feature ('SpatialPolygonsDataFrame' class) representing the
#' study area and exactly matching the extent of the covariate layer. If the latter is not provided,
#' it is internally worked out from the covariate raster and may make the whole function take a
#' while to complete.\cr
#'
#' The function plots a frequency distribution histogram of the average value of a spatial covariate at
#' randomized locations (using B iterations). At each iteration, the number of randomized points is
#' equal to the number of points of the input point pattern. Two blue reference lines correspond
#' to the 0.025th and to the 0.975th quantile of the randomized distribution.
#' A black dot represents the observed mean value of the covariate measured at the locations of
#' the input point pattern. P-values are reported. \cr
#'
#' @param feature Feature (of point type; 'SpatialPointsDataFrame' class) representing the spatial
#'   point pattern of interest.
#' @param cov.var Numeric covariate (of 'RasterLayer' class).
#' @param studyplot Feature (of polygon type; 'SpatialPolygonsDataFrame' class) representing the
#'   study area and exactly matching the extent of the covariate layer. If NULL, it is worked out
#'   from the covariate layer (may make the whole function take a while to complete).
#' @param B Number of randomized iterations to be used to calculate the acceptance interval (199 by
#'   default).
#' @param oneplot Set to TRUE (default), will plot the charts into a single visualization.
#'
#' @keywords pointsCovarDistr
#'
#'@return The function returns a list storing the following components: \itemize{
##'  \item{$obs.cov.values: }{observed values of the covariate at the point pattern locations}
##'  \item{$obs.average: }{average of the observed values of the covariate}
##'  \item{$p.value.obs.smaller.than.exp: }{p.value for the observed average smaller than expected under the Null Hypothesis}
##'  \item{$p.value.obs.larger.than.exp: }{p.value for the observed average larger than expected under the Null Hypothesis}
##'  \item{$p.value.obs.diff.from.exp: }{p.value for the observed average different from what expected under the Null Hypothesis}
##' }
#'
#' @export
#'
#' @importFrom raster rasterToPolygons extract
#'
#' @examples
#' #load the point dataset representing the location of springs
#' data(springs)
#'
#' #load the polygon dataset representing the study area
#' data(malta_polyg)
#'
#' #load the raster representing the terrain elevation, to be used as covariate
#' data(malta_dtm_40)
#'
#'pointsCovarDistr(feature=springs, cov.var=malta_dtm_40, studyplot=malta_polyg)
#'
#' @seealso \code{\link{pointsCovarModel}}, \code{\link{pointsCovarCum}}
#'
#'
pointsCovarDistr <- function (feature, cov.var, studyplot = NULL, B = 199, oneplot = TRUE){
  #disable scientific notation
  options(scipen = 999)

  #if the studyplot is not provided, create the polygon representing the studyplot using the
  #covariate raster
  if (is.null(studyplot) == TRUE) {
    cov.var.copy <- cov.var
    cov.var.copy[na.omit(cov.var.copy)] <- 1
    studyplot <- raster::rasterToPolygons(cov.var.copy, na.rm = TRUE, dissolve = TRUE)
  }
  else {
  }

  region <- studyplot

  #extract the value of the covariate at the points location
  cov.var.extract <- raster::extract(cov.var, feature)

  #calculate the aboserved mean of the covariate at the points location
  cov.var.extract.mean <- mean(cov.var.extract, na.rm=TRUE)

  #create an empty vector of length B+1
  mean.cov.rnd <- length(B+1)

  #store the observed mean of the covariate in the first slot of the empty vector
  #created in the preceding step
  mean.cov.rnd[1] <- cov.var.extract.mean

  #set up the progress bar
  pb <- txtProgressBar(min = 0, max = B, style = 3)

  #loop to draw a random sample of points within the study plot, calculating the mean of the covariate
  #at the random points location, and storing sequentialy in the vector created in one of the preceding
  #steps
  for (i in 2:B+1) {
    rnd <- spsample(region, n = length(feature), type = "random")
    cov.var.rnd <- raster::extract(cov.var, rnd)
    mean.cov.rnd[i] <- mean(cov.var.rnd, na.rm=TRUE)
    setTxtProgressBar(pb, i)
  }

  #calculate the p values
  p.lowertail <- (1 + sum (mean.cov.rnd[-1] < mean.cov.rnd[1], na.rm=TRUE)) / (1 + B)
  p.uppertail <- (1 + sum (mean.cov.rnd[-1] > mean.cov.rnd[1], na.rm=TRUE)) / (1 + B)
  two.sided.p <- 2 * min(p.lowertail, p.uppertail)


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

  #set the title for the histogram
  maintitle <- paste0("Frequency distribution of the average covariate value at random points' location \n(based on ", B, " randomized iterations)")

  #plot the frequency distribution histogram
  graphics::hist(mean.cov.rnd, main=maintitle,
                 sub=paste0("Observed mean: ", round(mean.cov.rnd[1], 3),"\n p-value obs. average smaller than expected: ", p.lowertail, "; p-value obs. average larger than expected: ", p.uppertail, "\np-value obs. average different from expected: ", two.sided.p),
                 xlab="",
                 cex.main=0.90,
                 cex.sub=0.70,
                 cex.axis=0.85)

  #add the rug plot
  rug(mean.cov.rnd, col = "#0000FF")

  #plot the vertical lined representing lower and upper tail
  abline(v=stats::quantile(mean.cov.rnd, 0.025, na.rm=TRUE), lty=2, col="blue")
  abline(v=stats::quantile(mean.cov.rnd, 0.975, na.rm=TRUE), lty=2, col="blue")

  #plot the abserved mean value of the covariate at the points' location
  points(x=mean.cov.rnd[1], y=0, pch=20, col = "black")

  #restore the original graphical device's settings
  par(mfrow = c(1,1))

  results <- list("obs.cov.values"=cov.var.extract,
                  "obs.average" = cov.var.extract.mean,
                  "p.value.obs.smaller.than.exp"=round(p.lowertail,3),
                  "p.value.obs.larger.than.exp"=round(p.uppertail,3),
                  "p.value.obs.diff.from.exp"=round(two.sided.p,3))

  return(results)
}
