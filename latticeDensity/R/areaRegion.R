#' Computes area of a region
#' 
#' This function computes the area of a region by first finding
#' the area of the bounding polygon, then subtracting the area
#' of each hole.
#' 
#' @section Warning:
#' Note that this program does not check to see if the holes are non-intersecting
#' or if the holes intersect the polygon.
#' @param formLatticeOutput An object returned by formLattice or
#' editLattice.
#' @examples 
#' data(areaRegionExample)
#' attach(areaRegionExample)
#' hole_list <- list(hole1,hole2)
#' nodeFillingOutput <- nodeFilling(poly=boundary, node_spacing=0.03,
#'                                hole_list = hole_list)
#' formLatticeOutput <- formLattice(nodeFillingOutput)
#' plot(formLatticeOutput)
#' areaRegion(formLatticeOutput)
#' @import utils
#' @import graphics
#' @import stats
#' @import splancs
#' @export
areaRegion <-
function(formLatticeOutput){
  #
  #  Computer the area of the bounding polygon
  #  and subtract the areas of the hole polygons
  #
boundary <- formLatticeOutput$poly
hole_list <- formLatticeOutput$hole_list
number_holes <- length(hole_list)
hole_areas <- rep(NA, number_holes)
if(number_holes > 0){
  for (i in 1:number_holes){
    hole_areas[i] <- splancs::areapl(hole_list[[i]])
  }
  areaRegion <- splancs::areapl(boundary) - sum(hole_areas)
}
if(number_holes == 0){
  areaRegion <- splancs::areapl(boundary)
  }
return(areaRegion)
}
