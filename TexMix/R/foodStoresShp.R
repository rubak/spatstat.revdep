#' @title Point layer of Stores selling Food in Dallas County, TX
#' @description Location of food stores in Dallas County, TX, in the longitude and
#'   latitude format (see \code{proj4string=CRS("+proj=longlat +ellps=WGS84")})..
#' @docType data
#' @name foodStoresShp
#' @author The data were compiled by Laila Al Aswad and Michael Tiefelsdorf <tiefelsdorf@@utdallas.edu>
#' @source Reference USA, 2019 \url{http://www.referenceusa.com}.
#' @format Spatial polygon data-frame with 1623 verified store locations.
#' \describe{
#'   \item{SALESVOL}{Reported total annual sales volume of goods in $}
#'   \item{PROPFOOD}{Assumed proportion of food sales}
#'   \item{FOODSALES}{Calculated annual sales volume of food in $}
#'   \item{STORETYPE}{Factor distinguishing between stores selling nutritious
#'                    food (grocery stores) and processed food (convenience
#'                    stores)}
#' }
#' @examples
#' library(spatstat.core)
#' library(rgdal)
#' library(sp)
#' proj4string(bndShp)                                     # Current system
#' projUTM <- CRS("+proj=utm +zone=14  +units=m")          # isotropic coordinate sytem
#' bndUTM <- spTransform(bndShp, projUTM)                  # Re-project boundary
#' storesUTM <- spTransform(foodStoresShp, projUTM)        # Re-project points
#' storesDf <- as.data.frame(storesUTM)                    # Extract data-frame
#' storesPts <- as.ppp(storesUTM)                          # Convert to .ppp
#' storesPts$marks <- NULL                                 # Clear marks
#' bndWin <- as.mask(as.owin(bndUTM), eps=200)             # pixel window with 200 m resolution
#' unitname(bndWin) <- list("meter","meters")              # set units
#' storesPts <- storesPts[bndWin]                          # assign window to pts
#' summary(storesPts)
#'
#' ## Evaluate weighted kernel density with bw=3000
#' allFoodIm <- density(storesPts, weights=storesDf$FOODSALES, sigma=3000)
#' plot(allFoodIm, main="All Stores Weighted Kernel Density\nbw = 3000 m")
#' plot(storesPts, cex=0.5, pch=16, col="green", add=TRUE)
#' box(); axis(1); axis(2)
#'
NULL
