#' R function to model (and test) the dependence of a point pattern on a spatial numeric covariate
#'
#' The function is a wrapper for a number of functions out of the extremely useful 'spatstat'
#' package (specifically, ppm(), cdf.test(), auc(), roc(), effectfun()). It allows to test if there
#' is a significant dependence of the input point pattern on a underlying spatial numeric covariate
#' (first-order effect).\cr The function takes as input three datasets: a point patter
#' ('SpatialPointsDataFrame' class), a covariate layer (of 'RasterLayer' class), and a polygon
#' feature ('SpatialPolygonsDataFrame' class) representing the study area and exactly matching the
#' extent of the covariate layer. If the latter is not provided, it is internally worked out from
#' the covariate raster and may make the whole function take a while to complete.\cr
#'
#' The function fits a inhomogeneous Poisson point process (Alternative Model-H1) with intensity of
#' the point pattern as a loglinear function of the underlaying numerical covariate (see Baddeley et
#' al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016, 307-309).
#' Also, the function fits a homogeneous Poisson point model (Null Model-H0, equivalent to Complete
#' Spatial Randomness: Baddeley et al., "Spatial Point Patterns. Methodology and Applications with
#' R", CRC Press 2016, 305-306), that is used as comparison for the inhomogeneous point process
#' model in a Likelihood Ratio test (Baddeley et al., "Spatial Point Patterns. Methodology and
#' Applications with R", CRC Press 2016, 334-335). A significant result, i.e. a low p-value,
#' suggests rejecting the Null Hypothesis of CSR in favour of the Alternative Hypothesis of a
#' Poisson point process affected by a covariate effect (i.e., inhomogeneous intensity due to the
#' influence of the covariate) (Baddeley et al., "Spatial Point Patterns. Methodology and
#' Applications with R", CRC Press 2016, 305). \cr
#'
#' @param feature Feature (of point type; SpatialPointsDataFrame class) representing the spatial
#'   point pattern of interest.
#' @param cov.var Numeric covariate (of RasterLayer class).
#' @param studyplot Feature (of polygon type; SpatialPolygonsDataFrame) representing the study area
#'   and exactly matching the extent of the covariate layer. If NULL, it is worked out from the
#'   covariate layer (may make the whole function take a while to complete).
#' @param oneplot Set to TRUE (default), will plot the charts into a single visualization.
#'
#' @return  The function returns a 4 plots, which can be arranged in just one visualization setting
#' the parameter 'oneplot' to TRUE:\cr
#'
#' -plot of the point pattern along with the underlaying covariate raster;\cr
#'
#' -plot of the fitted intensity against the spatial covariate (Baddeley et al., "Spatial Point
#' Patterns. Methodology and Applications with R", CRC Press 2016, 308);\cr
#'
#' -plot of the cumulative distribution of the covariate at the data points against the cumulative
#' distribution of the covariate at all the spatial location within the study area (rationale:
#' Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016,
#' 184-185);\cr
#'
#' -plot of the ROC curve, which help assessing the strength of the dependence on the covariate
#' (Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016,
#' 187-188).\cr
#'
#' -a list containing: \itemize{
##'  \item{$H0-model: }{info and relevant statistics regarding the Null Model}
##'  \item{$H1-model: }{info and relevant statistics regarding the Alternative Model}
##'  \item{$Model comparison (LRT): }{results of the Likelihood Ratio test}
##'  \item{$AIC-H0: }{AIC of the Null Model}
##'  \item{$AIC-H1: }{AIC of the Alternative Model}
##'  \item{$KS test: }{information regarding the cumulative distribution comparison via
##'  Kolmogorov-Smirnov test}
##'  \item{$AUC: }{AUC statistic}
##' }
#'
#' @keywords covariate
#'
#' @export
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
#' #note: a warning message reporting that 4 out of 699 points
#' # have values of the covariate undefined is expected
#' results <- pointsCovarModel(Starbucks, popdensity, Massachusetts)
#'
#' @seealso \code{\link{pointsCovarCum}} , \code{\link{distCovarModel}} , \code{\link{distRandSign}}
#'
pointsCovarModel <- function(feature, cov.var, studyplot=NULL, oneplot=FALSE){
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

  #transform the feature variable to a ppp class (needed by spatstat) and remove appended data by using ummark()
  feature.ppp <- unmark(as.ppp(feature))

  #set the analytical window (needed by spatstat) to the extent of the studyplot
  W  <- as(studyplot, "owin")

  #give the feature dataset the same analytical window of the study region
  spatstat::Window(feature.ppp) <- W

  #tranform the cov.var from a RasterLayer to an object of class im, which is needed by spatstat
  cov.var.im <- spatstat::as.im(cov.var)

  #fit the Null Model (homegeneous point process)
  PPM0 <- ppm(feature.ppp ~ 1)

  #fit the Alternative Model (inhomogeneous point process with intensity as function of the cavariate)
  PPM1 <- ppm(feature.ppp ~ cov.var.im)

  #perform the KS test
  kolmsmirn <- cdf.test(feature.ppp, cov.var.im)

  #calculate the AUC
  areaundercurve <- spatstat::auc(feature.ppp, cov.var.im, high=FALSE)

  #compare the models via likelihood ratio test
  model.comp <- anova(PPM0, PPM1, test="LRT")

  #set the output of the graphic device according to the oneplot parameter
  if(oneplot==TRUE){
    par(mfrow=c(2,2))
  } else {}

  #extract the p-value of the likelihod ratio test to be used in the plot subtitle
  anova.p <- model.comp$"Pr(>Chi)"[2]

  #classify the anova p value to be used in the plot subtitle
  anova.p.to.report <- ifelse(anova.p < 0.001, "< 0.001",
                              ifelse(anova.p < 0.01, "< 0.01",
                                     ifelse(anova.p < 0.05, "< 0.05",
                                            round(anova.p, 3))))

  #plot the study region along with the point pattern and the feature whose distance is used as covariate
  raster::plot(cov.var, main="Point pattern against the numeric covariate data", cex.main=0.75,
       sub=paste0("Null Hypothesis (H0): Homogeneous Poisson process model\nAlternative Hyphotesis (H1): Inhomogeneous Poisson process model (intensity as loglinear function of the covariate)\nH1 ", ifelse(anova.p > 0.05, "is not", "is"), " a significant improvement over H0 (Likelihood Ratio p-value: ", anova.p.to.report,"; AUC: ", round(areaundercurve,3), ")"), cex.sub=0.70)

  raster::plot(feature, add=TRUE, pch=20)
  raster::plot(studyplot, add=TRUE)

  #plot the fitted Alternative Model
  #i.e., modelled intensity against the covariate
  plot(spatstat::effectfun(PPM1, names(PPM1$covariates), se.fit=TRUE),
       main="Fitted intensity of the point pattern \nas (loglinear) function of the covariate",
       cex.main=0.8,
       cex.axis=0.7,
       cex.lab=0.8, legend=TRUE)

  #plot the cumulative distr chart
  plot(kolmsmirn, cex.main=0.8)

  #plot the ROC curve
  plot(spatstat::roc(PPM1),
       main=paste0("ROC curve of the fitted intensity of point patter \nas (loglinear) function of the cavariate \nAUC: ", round(areaundercurve,3)),
       cex.main=0.8)

  #create a list to store relevant results
  results <- list("H0-model"=PPM0,
                  "H1-model"=PPM1,
                  "Model comparison (LRT)"=model.comp,
                  "AIC-H0"=AIC(PPM0),
                  "AIC-H1"=AIC(PPM1),
                  "KS test"=kolmsmirn,
                  "AUC"=areaundercurve)

  # restore the original graphical device's settings
  par(mfrow = c(1,1))

  return(results)
}
