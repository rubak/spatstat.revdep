#' @title Terrain Ruggedness Index
#' @description Implementation of the Riley et al (1999) Terrain Ruggedness Index
#'
#' @param r              RasterLayer class object
#' @param s              Scale of window. Must be odd number, can represent 2 
#'                       dimensions (eg., s=c(3,5) would represent a 3 x 5 window)
#' @param exact          Calculate (TRUE/FALSE) the exact TRI or an algebraic 
#'                       approximation. 
#' @param file.name      Name of output raster (optional)
#' @param ...            Additional arguments passed to writeRaster
#'
#' @return raster class object or raster written to disk
#'
#' @description
#' The algebraic approximation is considerably faster. However, because 
#' inclusion of the center cell, the larger the scale the larger the divergence 
#' of the minimum value.
#' 
#' @description
#' Recommended ranges for classifying Topographic Ruggedness Index:
#' * 0-80 - level terrain surface.
#' * 81-116 - nearly level surface.
#' * 117-161 - slightly rugged surface.
#' * 162-239 - intermediately rugged surface.
#' * 240-497 - moderately rugged surface.
#' * 498-958 - highly rugged surface.
#' * gt 959 - extremely rugged surface. 
#' @md 
#'
#' @author Jeffrey S. Evans  <jeffrey_evans@@tnc.org>
#'  
#' @references 
#' Riley, S.J., S.D. DeGloria and R. Elliot (1999) A terrain 
#'   ruggedness index that quantifies topographic heterogeneity, 
#'   Intermountain Journal of Sciences 5(1-4):23-27.
#'
#' @examples 
#' \donttest{
#'  library(raster)
#'  data(elev)
#'   ( tri.ext <- tri(elev) )
#'   ( tri.app <- tri(elev, exact = FALSE) )
#'   plot(stack(tri.ext, tri.app))
#' }
#'
#' @export
tri <- function(r, s = 3, exact = TRUE, file.name = NULL, ...) {
  if(length(s) > 2) stop( "Specified window exceeds 2 dimensions")   
    if(length(s) == 1) s = rep(s,2)
  r.sqrt <- function(x) {
    x <- x[x >= 0]
    if(length(x) >= 1) { 
      return(sqrt(x)) 
    } else{
      return(NA)  
    }
  }	
  tri.calc <- function(x) {
    xc <- x[ceiling(length(x)/2)]
    x <- x[-ceiling(length(x)/2)]
    x <- x[!is.na(x)]
      if(!is.na(xc) & length(x) > 0) {
        x.dev <- vector()
        for(i in 1:length(x)) { x.dev <- append(x.dev, (xc - x[i])^2) }
  	  return( sqrt(sum(x.dev)) )
      } else {
      return( NA )  
    }
  }	  
  if(exact == TRUE) {  
    if(!is.null(file.name)) {
      return( raster::focal(r, w=matrix(1,s[1],s[2]), fun=tri.calc, na.rm=FALSE,
	                        filename = file.name, ...) )
    } else {
      return( raster::focal(r, w=matrix(1,s[1],s[2]), fun=tri.calc, na.rm=FALSE) )
    }  
  } else {
    sx <- raster::focal(r, w = matrix(1,s[1],s[2]), sum)
    e2 <- r * r
    st <- raster::focal(e2, w = matrix(1,s[1],s[2]), sum)
    r2 <- st + (s[1]*s[2]) * e2 - 2 * r * sx
	if(!is.null(file.name)) {
      return( raster::calc(r2, fun= r.sqrt, filename = file.name, ...) )
	} else {
	  return( raster::calc(r2, fun= r.sqrt) )
    }
  }
}
