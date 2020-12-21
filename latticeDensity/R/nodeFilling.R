#' Produces a grid of locations inside the region boundary.
#' 
#' nodeFilling produces a grid of locations that are the nodes 
#' in the diffusion process.
#' 
#' @details nodeFilling superimposes a square grid 
#' of points over the region, with spacing given 
#' by the parameter node_spacing. The points 
#' contained in the region are retained. The output, 
#' a nodeFillingOutput object, contains the boundaries 
#' of the region (and holes), the set of nodes, and EW
#' and NS coordinates necessary for creating a contour plot.
#' 
#' @param poly A matrix that contains the vertices of the 
#' bounding polygon.
#' @param node_spacing The distance between grid locations.
#' @param hole_list Optional list of holes to be removed from 
#' the region
#' @return An object of type nodeFillingOutput is produced.
#' \itemize{
#'   \item EW_locs EW coordinates for the contour plot.
#'   \item NS_locs NS coordinates for the contour plot.
#'   \item nodes Matrix of node locations.
#'   \item poly Matrix of vertices of boundary polygon.
#'   \item node_spacing Vertical and horizontal node spacing.
#'   \item hole_list List of polygons representing holes in region.}
#' @author Ronald P. Barry
#' @references Ronald P. Barry, Julie McIntyre. Estimating 
#' animal densities and home range in regions with irregular 
#' boundaries and holes: A lattice-based alternative to the 
#' kernel density estimator. Ecological Modelling 222 (2011) 
#' 1666-1672.
#' @examples 
#' plot.new()
#' data(polygon1)
#' nodeFillingOutput <- nodeFilling(poly=polygon1,node_spacing=0.02)
#' plot(nodeFillingOutput)
#' formLatticeOutput <- formLattice(nodeFillingOutput)
#' plot(formLatticeOutput)
#' 
#' @import utils
#' @import graphics
#' @import stats
#' @import splancs
#' @export
nodeFilling <- function(poly, node_spacing,hole_list = NULL){
   poly <- as.matrix(poly)
   if(!is.null(hole_list)){
     number.holes <- length(hole_list)
     for(k in 1:number.holes){
       hole_list[[k]] <- as.matrix(hole_list[[k]])}
     print("This function does not check to see if the holes")
     print("are nonintersecting, or whether they are contained")
     print("inside the boundary")
   }
   node_spacing <- as.numeric(node_spacing)
   #
   #  Fill the region with a grid of nodes
   #
   EW_locs <- seq(min(poly[,1])+node_spacing/2,
                  max(poly[,1])-node_spacing/2,by=node_spacing)
   NS_locs <- seq(min(poly[,2])+node_spacing/2,
                  max(poly[,2])-node_spacing/2,by=node_spacing)
   bound_array <- expand.grid(EW_locs,NS_locs)
   names(bound_array) <- c("x","y")
   nodes <- bound_array[splancs::inout(pts=bound_array, 
                        poly=poly, bound=TRUE),]
   nodes <- splancs::as.points(as.matrix(nodes))
   #
   #
   #
   output <- list(EW_locs = EW_locs,
                 NS_locs = NS_locs,
                 nodes = nodes,
                 poly = poly,
                 node_spacing = node_spacing,
                 hole_list = hole_list)
   class(output) <- "nodeFillingOutput"
   if(!is.null(hole_list)){for (k in 1:number.holes){
     output <- removeHole(hole_list[[k]], output)}}
   return(output)
}


