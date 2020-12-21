
#' @importFrom graphics   arrows axis par points text
#' @importFrom methods    as new slot
#' @importFrom stats      rnorm rpois

#' @importFrom rgeos      gBuffer gIntersection gUnionCascaded gTouches readWKT
#' @importFrom spdep      poly2nb card nblag
#' @importFrom deldir     deldir tile.list
#' @importFrom spatstat   ripras

#' @import sp   



### 

.onAttach <- function(libname, pkgname) {
	dcf <- read.dcf(file=system.file("DESCRIPTION", package=pkgname) )
    
	packageStartupMessage(paste('This is', pkgname, dcf[, "Version"], 'For a tutorial type vignette("expp")' ))
  

	
}

