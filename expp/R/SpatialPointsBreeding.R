#' @rdname      SpatialPointsBreeding
#' @exportClass  SpatialPointsBreeding
setClass("SpatialPointsBreeding", representation(
    id          = "numeric",
    male        = "character", 
    female      = "character"
    ),
    
    contains  = "SpatialPointsDataFrame",

    validity = function(object) {

        if (length(table(object@id)[table(object@id) > 1]) )
            stop("only one id per line is allowed.")
        if( any(object@id < 1) )            
            stop("id < 1 not allowed.")
        
        return(TRUE)
        # TODO
            # polygynous males
            # multiple br. att.
        }     
 )


#' Converts a \code{data.frame} to an object of class
#' \code{SpatialPointsBreeding}
#' 
#' Converts a \code{\link{data.frame}} to a \code{SpatialPointsBreeding}
#' object. The \code{SpatialPointsBreeding} class extends
#' \code{\link[sp]{SpatialPointsDataFrame}} with three extra slots defining the
#' id (i.e. nest or breeding box) and the pair identity (i.e. male and female),
#' respectively.
#' 
#' @param data a \code{\link{data.frame}} containing the coordinates 
#'        (e.g. "x","y"), the location id, and the pair identity 
#'        (e.g. "male", "female") together with any other optional variables 
#'        (e.g. individuals or nest traits).
#' @param proj4string A \code{\link[sp]{CRS}} object containing a valid proj4
#'        string.  See \code{\link[sp]{CRS}} \code{\link[sp]{proj4string}} 
#'        for details.
#' @param coords Formula specifying which columns in object are the spatial
#'        coordinates.  Argument passed to \code{\link[sp]{coordinates}}
#' @param breeding One side formula defining the male and female ID in that order (e.g. ~ male + female)
#' @param id Integer specifying the location id (e.g. nest box number, den ID).
#' @param x   a \code{SpatialPointsBreeding} object
#' @param y   an \code{eppMatrix} object
#' @param pch see \code{plot.default}
#' @param axes see \code{plot.default}
#' @param add  see \code{plot.default}
#' @param xlim see \code{plot.default}
#' @param ylim see \code{plot.default}
#' @param \dots further arguments to pass to plot(as(x, "Spatial")
#' @param cex see \code{plot.default}
#' @param col see \code{plot.default}
#' @param col.epp extra-pair partners color
#' @param lwd see \code{plot.default}
#' @param lty see \code{plot.default}
#' @param bg see \code{plot.default}
#' @export
#' @include eppMatrix.R
#' @seealso \code{\link[expp]{epp}}
#' @examples
#' d = data.frame(
#'   x = c(4, 17, 16, 41, 41, 43, 86, 62, 71, 92, 95,53, 34, 27, 53), 
#'   y = c(3, 18, 36, 6, 18, 50, 3, 21, 40, 43, 57, 62, 62, 45, 37), 
#'   id = 1:15,male = paste0('m', 1:15), female = paste0('f', 1:15), 
#'     stringsAsFactors = FALSE)
#' 
#' b = SpatialPointsBreeding(d, id = 'id', breeding = ~ male+female)    
#' 
#' plot(b)
#' 
#' 
SpatialPointsBreeding <- function(data,  proj4string, coords = ~ x + y, breeding = ~ male + female, id ) {
	d = data
	row.names(d) = NULL
	d$k = 1:nrow(d)
	coordinates(d) <- coords
    if(missing(proj4string)) proj4string  = CRS(as.character(NA))
	   proj4string(d) = proj4string
	
	ids = data[, id]
	
	m = as.character(breeding[[2]][2])
	f  = as.character(breeding[[2]][3])
	males = data[, m]
	females = data[, f]
	
	d@data[, m]  = NULL
	d@data[, f]  = NULL
	d@data[, id] = NULL
	
	new("SpatialPointsBreeding", d, id = ids, male = males, female= females)
    }



if (!isGeneric("plot")) setGeneric("plot", function(x, y, ...) standardGeneric("plot"))

#' @rdname 	SpatialPointsBreeding
setMethod("plot", signature(x = "SpatialPointsBreeding", y = "missing"),
          function(x, pch = 20, axes = FALSE, add = FALSE, 
                   xlim = NULL, ylim = NULL, ..., cex = 1, col = "grey", lwd = 1, bg = "grey90") {
            if (! add)
              plot(as(x, "Spatial"), axes = axes, xlim = xlim, ylim = ylim, ...)
            cc = coordinates(x)
            points(cc[,1], cc[,2], pch = pch, cex = cex, col = col, lwd = lwd, bg = bg)
            text(cc[,1], cc[,2], x@id, pos = 4, cex = cex)
            text(cc[,1], cc[,2], x@female, pos = 1,  cex =  cex-0.1)
            text(cc[,1], cc[,2], x@male, pos = 3,  cex = cex-0.1)
            
          })

#' @rdname  SpatialPointsBreeding		  
setMethod("plot", signature(x = "SpatialPointsBreeding", y = "eppMatrix"),
          function(x, y, pch = 20, axes = FALSE, add = FALSE, 
                   xlim = NULL, ylim = NULL, ..., cex = 1, col = "grey", col.epp = "red", lwd = 1, lty = 2, 
                   bg = "grey90") {
            if (! add)
              plot(as(x, "Spatial"), axes = axes, xlim = xlim, ylim = ylim, ...)
            cc = coordinates(x)
			
			# nests
			points(cc[,1], cc[,2], pch = pch, cex = cex, col = col, lwd = lwd, bg = bg)
            text(cc[,1], cc[,2], x@id, pos = 4, cex = cex)
            
			# males
			epm = which(x@male %in% y@male)
			try(text(cc[epm,1], cc[epm,2], x@male[epm], pos = 1,  cex =  cex-0.1, col = col.epp), silent = TRUE)
			try(text(cc[-epm,1], cc[-epm,2], x@male[-epm], pos = 1,  cex =  cex-0.1), silent = TRUE)
			
			# females
			epf = which(x@female %in% y@female)
			try(text(cc[epf,1], cc[epf,2], x@female[epf], pos = 3,  cex =  cex-0.1, col = col.epp), silent = TRUE)
			try(text(cc[-epf,1], cc[-epf,2], x@female[-epf], pos = 3,  cex =  cex-0.1), silent = TRUE)
			
			# connections
			for(i in 1:length(y@male) ) {
				mc = cc[which(x@male == y@male[i]), ]
				if(length(mc) == 0) warning("EP male ", sQuote(y@male[i]), " not found.")
				if( length(mc) > 2) mc = mc[1, ] # only one line per polygynous male
					
				fc = cc[which(x@female == y@female[i]), ]
				if(length(fc) == 0) warning("EP female ", sQuote(y@female[i]), " not found.")
								
				arrows(mc[1], mc[2], fc[1], fc[2], col = col.epp, code = 3, angle = 12, length = 0.2, lty = lty)
				}
           
          })
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  

	











