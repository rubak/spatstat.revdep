#' @rdname      eppMatrix
#' @slot    male   extra-pair male   ID-s as character vectors
#' @slot    female extra-pair female ID-s as character vectors
#' @exportClass eppMatrix
setClass("eppMatrix", representation(
    male   = "character", 
    female  = "character" 
    ),

    validity = function(object) {
    if( length( intersect(object@male, object@female) ) > 0 ) stop("the same id cannot be male and female in the same time .")
    if ( any( is.na(object@male) ) ) stop("NA values are not allowed.")
    if ( any( is.na(object@female) ) ) stop("NA values are not allowed.")

    return(TRUE)
    }
    )


#' Convert a \code{data.frame} to an eppMatrix object.
#' 
#' Converts a \code{data.frame} to a eppMatrix object using a
#' \code{~male+female} formula.
#' 
#' @param data a \code{data.frame}
#' @param pairs a formula indicating the extra-pair male and the extra-pair female in that order.
#' @return An object of class \code{eppMatrix} with two slots.

#' @seealso \code{\link[expp]{epp}}
#' @export 
#' @examples
#' 
#' eppPairs = data.frame(male = c("m1", "m2", "m1"), female=c("f3", "f1", "f2") )
#' e = eppMatrix(eppPairs,  pairs = ~ male + female)
#' class(e)
#' showClass("eppMatrix")
#' 
#' data(bluetit_breeding)
#' data(bluetit_epp)
#' b = bluetit_breeding[bluetit_breeding$year_ == 2010, ]
#' eppPairs = bluetit_epp[bluetit_epp$year_ == 2010, ]
#' 
#' breedingDat  = SpatialPointsBreeding(b, id = 'id', coords = ~ x + y, breeding = ~ male + female)
#' eppDat = eppMatrix(eppPairs, pairs = ~ male + female)
#' 
#' plot(breedingDat, eppDat)
#' 
eppMatrix <- function(data,  pairs = ~ male + female) {
	
	m = as.character(pairs[[2]][2])
	f  = as.character(pairs[[2]][3])
  
  # TODO: remove repeated lines, warn
  
  new('eppMatrix', male = as.character(data[, m]), female = as.character(data[, f]) )
	
  }













