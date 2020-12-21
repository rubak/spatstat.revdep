#' Dirichlet Polygons
#' 
#' Computes the Dirichlet polygons using a
#' \code{\link[expp]{SpatialPointsBreeding}} object and optionally a boundary
#' \code{\link[sp]{SpatialPolygons}} or a vector containing id-s located at the boundary.
#' @param x         A \code{\link{SpatialPointsBreeding}} object.
#' @param boundary  A \code{\link[sp]{SpatialPolygons}} or a vector of integers containing the id-s located at the boundary. 
#'                  When missing boundary is inferred using \code{\link[spatstat.geom]{ripras}} in \code{spatstat} .
#' @param width     argument passed to \code{\link[rgeos]{gBuffer}}. It defines the distance between boundary boxes and the boundary polygon; 
#'                  it is set by default to half of the average distance between boundary boxes.
#' @param \dots     passed to \code{\link[spatstat.geom]{ripras}}
#' @export
#' @include SpatialPointsBreeding.R
#' @examples
#' 
#' d = data.frame(
#'   x = c(4, 17, 16, 41, 41, 43, 86, 62, 71, 92, 95,53, 34, 27, 53), 
#'   y = c(3, 18, 36, 6, 18, 50, 3, 21, 40, 43, 57, 62, 62, 45, 37), 
#'   id = 1:15,male = paste0('m', 1:15), female = paste0('f', 1:15), 
#'     stringsAsFactors=FALSE)
#' 
#' b = SpatialPointsBreeding(d, id = 'id', breeding = ~ male + female)  
#' 
#' # boundary is inferred based on the Ripley-Rasson estimate of the spatial domain
#' dp1 = DirichletPolygons(b)
#' plot(dp1)
#' 
#' # boundary is given
#' brdy2 = rgeos::readWKT("POLYGON((28 71,67 68,70 49,84 49,90 74,111 65,107 
#'                       36,78 28,98 15,98 -4,74 -7,-2 -8,0 31,28 71) )")
#' dp2 = DirichletPolygons(b, boundary = brdy2)
#' plot(dp2)
#' 
#' # boundary is inferred based on the boundary id-s. 
#' # define boundary id-s using a 'Follow-The-Dots' strategy. 
#' brdy3 = as.integer(c(1, 2, 4, 7, 9, 10, 11, 12, 13, 14, 3))
#' 
#' dp3 = DirichletPolygons(b, boundary = brdy3)
#' plot(dp3)
#' 
#' # setting width manually
#' dp4 = DirichletPolygons(b, boundary = brdy3, width = 2)
#' plot(dp4)
#' 
#' 
#' plot(dp1)
#' plot(dp2, add = TRUE, border = 2)
#' plot(dp3, add = TRUE, border = 3)
#' plot(dp4, add = TRUE, border = 4)
#' plot(b, add = TRUE)
#' 
#' 


#' @rdname      DirichletPolygons
setGeneric("DirichletPolygons", function(x, boundary, ...)   standardGeneric("DirichletPolygons") )


.DirichletPolygons <- function(x, boundary) {
		# x is a SpatialPointsBreeding, boundary is a SpatialPolygons*
      coords = coordinates(x)
            
			p  =  tile.list(deldir(coords[,1], coords[,2], suppressMsge = TRUE))
			p  = lapply(p, function(P) data.frame(x = P$x, y = P$y) )
			
			d  =  do.call(rbind, p)    
			d$id = rep(x@id, sapply(p, nrow) )
			
			P = split(d, d$id)
			P = lapply(P, function(x) { x = 
							x = rbind(x, x[1, ])  
							row.names(x) = paste(x$id, 1:nrow(x), sep = "_")
							x
							} )
			
			P = lapply(P,  function(x) SpatialPolygons(list( Polygons(list(Polygon(x[, c('x', 'y')])), ID = x$id[1] ) )) )

			P = lapply(P, function(pp) { proj4string(pp) = CRS(proj4string(x)); pp } )
      
			P = lapply(P,function(pj) gIntersection(boundary, pj, id = slot(slot(pj, 'polygons')[[1]], 'ID')))

			P = do.call(rbind, P) 
			P = SpatialPolygonsDataFrame(P, data = data.frame(ID = x@id, row.names = x@id))
			
      # bdry line
      bl = as( gUnionCascaded(P), 'SpatialLines')

      P$outerPoly = as.vector(gTouches( P, bl , byid = TRUE))

      P
    }

#' @rdname      DirichletPolygons
#' @export
setMethod("DirichletPolygons",  
          signature  = c(x = "SpatialPointsBreeding", boundary = "missing"), 
          definition = function(x, ...) {
            coords = coordinates(x)
            ids = x@id
            rr = spatstat.geom::ripras(coords, shape = "convex", ...)
            rr = cbind(x = rr$bdry[[1]]$x, y = rr$bdry[[1]]$y)
            boundary =  SpatialPolygons(list( Polygons(list( Polygon(rbind(rr, rr[1, ] )) ) , 1) ) )
            proj4string(boundary) = proj4string(x)
            
            .DirichletPolygons(x, boundary)
          })

#' @rdname      DirichletPolygons
#' @export
setMethod("DirichletPolygons",  
          signature  = c(x = "SpatialPointsBreeding", boundary = "integer"), 
          definition = function(x, boundary, width) {
            
                z = data.frame(coordinates(x), id = x@id )

                bb = data.frame(id = boundary, o = 1:length(boundary) )
                bb = merge(bb, z, by = 'id')
                bb = bb[order(bb$o), ]
                bb = rbind(bb, bb[1, ])

                P = readWKT( paste( "POLYGON((", paste(paste(bb$x, bb$y), collapse = ','), "))" ) )

                if( missing(width) ) {
                        # median distance between points
                        z12 =  cbind(bb[-nrow(bb), c('x', 'y')], bb[-1, c('x', 'y')])
                        width = mean(apply(z12, 1, function(x) spDists( as.matrix(t(x[1:2])), as.matrix(t(x[3:4])) ) ) )/2
                        }
                
                P = gBuffer(P, width = width)         
                
      				.DirichletPolygons(x, P)
              })

#' @rdname      DirichletPolygons
#' @export
setMethod("DirichletPolygons",  
    signature  = c(x = "SpatialPointsBreeding", boundary = "SpatialPolygons"), 
    definition = function(x, boundary) {
		.DirichletPolygons(x, boundary)
    })


























