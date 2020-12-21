#' @title Lakes in Dallas County, TX
#' @description Polygon layer in the longitude and latitude format (see
#'   \code{proj4string=CRS("+proj=longlat +ellps=WGS84")}).
#' @docType data
#' @name lakesShp
#' @examples
#' library(maptools)
#' validTractShp <- tractShp[!is.na(tractShp$BUYPOW), ]  # Remove 2 tracts with NA's
#' plot(tractShp, col="white", border="white", axes=TRUE,
#'      main="Dallas Census Tracts with Food Deserts")
#' plot(validTractShp, col="ivory2", border="white", add=TRUE)
#' plot(lakesShp, col="skyblue", border="skyblue",add=TRUE)
#' plot(hwyShp, col="cornsilk3", lwd=3, add=TRUE)
#' plot(foodDesertShp, border="magenta",lwd=2, add=TRUE)
#' plot(bndShp, border="black", add=TRUE)
#' box()
#'
NULL
