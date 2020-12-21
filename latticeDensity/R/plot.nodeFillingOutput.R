#' Plot a nodeFillingOutput object.
#' 
#' Plots the boundary, all holes and 
#' the locations of all nodes. Should 
#' be used to decide if the nodes fill 
#' the region and are spaced closely 
#' enough to give good resolution in 
#' the plots. The only reason not to 
#' make the nodes too closely spaced 
#' is when the computing time or memory 
#' becomes too great.
#' 
#' @param x An object of type nodeFillingOutput 
#' returned by either nodeFilling or removeHole.
#' @param \dots Other arguments to be passed to 
#' functions plot, points, lines.
#' @author Ronald P. Barry
#' @references Ronald P. Barry, Julie McIntyre. 
#' Estimating animal densities and home range in 
#' regions with irregular boundaries and holes: 
#' A lattice-based alternative to the kernel density 
#' estimator. Ecological Modelling 222 (2011) 1666-1672.
#' @examples 
#' plot.new()
#' data(polygon1)
#' nodeFillingOutput <- nodeFilling(poly=polygon1, node_spacing=0.01)
#' plot(nodeFillingOutput)
#' 
#' @export
#' @method plot nodeFillingOutput
plot.nodeFillingOutput <-
function(x,...){
  if(class(x)!="nodeFillingOutput"){
       stop("Should be the output from the function nodeFilling")
    }
  nodes <- x$nodes
  poly <- x$poly
  hole_list <- x$hole_list
  wrapped_polygon <- rbind(poly,poly[1,])
  colnames(wrapped_polygon) <- c("Easting","Northing")
  plot(wrapped_polygon,lwd=2,type="l",...)
  points(nodes,pch=19,cex=0.5,...)
  if(!is.null(hole_list)){
    number_holes <- length(hole_list)
    for(k in 1: number_holes){
      lines(rbind(hole_list[[k]],hole_list[[k]][1,]),lwd=2,...)
    }
  }
}
