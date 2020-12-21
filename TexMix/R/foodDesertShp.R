#' @title Boundaries of three Neighorhoods in the City of Dallas, TX
#' @description Three polygons of putative food deserts in the longitude and
#' latitude format (see \code{proj4string=CRS("+proj=longlat +ellps=WGS84")}).
#' @docType data
#' @name foodDesertShp
#' @source Dallas Morning News Nov 16, 2016 (see
#' \url{https://www.dallasnews.com/opinion/editorials/2016/11/16/will-southern-dallas-food-deserts-get-relief})
#'
#' @format Spatial polygon data-frame with 3 areas and the following variables:
#' \describe{
#'   \item{ID}{Internal polygon ID}
#'   \item{AREA}{Neighborhood area in square miles}
#'   \item{DESERT}{Factor witht he name of the three neighborhoods}
#' }
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
