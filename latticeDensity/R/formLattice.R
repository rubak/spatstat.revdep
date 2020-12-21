#' Builds a neighbor structure on the nodes.
#' 
#' formLattice connects all nodes into a 
#' neighbor lattice by linking any two 
#' nodes that are within 1.5*node_spacing. 
#' Typically this will result in links in 
#' the E, W, N, S, NE, NW, SE, SW directions. 
#' The lattice object is created by the 
#' function dnearneigh from spdep.
#' 
#' @details 
#' When forming the lattice, the function does not 
#' check to see if any node is completely isolated 
#' from the rest of the nodes, nor does it check 
#' to see that paths exist between all pairs of nodes. 
#' Thus the lattice might be disconnected. You can 
#' still determine a nonparametric density in this 
#' case, but you need to think about whether it makes 
#' sense to allow disconnected sublattices. If you wish 
#' to connect isolated nodes to the lattice, use the 
#' editing function editLattice.
#' 
#' @param nodeFillingOutput An object, as produced by the function nodeFilling.
#' @return formLatticeOutput object
#' \itemize{
#'   \item EW_locs EW coordinates for use by contour
#'   \item NS_locs NS coordinates for use by contour
#'   \item nodes Matrix of node locations.
#'   \item poly Outer boundary.
#'   \item latt Neighbor lattice.
#'   \item hole.poly List of hole polygons.}
#' @author Ronald P. Barry
#' @references Ronald P. Barry, Julie McIntyre. 
#' Estimating animal densities and home 
#' range in regions with irreg- ular boundaries 
#' and holes: A lattice-based alternative to the 
#' kernel density estimator. Ecological Modelling 
#' 222 (2011) 1666-1672.
#' @examples 
#' plot.new()
#' data(polygon1)
#' nodeFillingOutput <- nodeFilling(poly=polygon1, node_spacing=0.02)
#' plot(nodeFillingOutput)
#' formLatticeOutput <- formLattice(nodeFillingOutput)
#' plot(formLatticeOutput)
#' Pointdata <- splancs::csr(polygon1,80)
#' Pointdata <- Pointdata[Pointdata[,1]<0.5,]
#' plot(polygon1,type="n")
#' polygon(polygon1)
#' points(Pointdata,pch=19)
#' densityOut <- createDensity(formLatticeOutput,PointPattern=Pointdata,
#'                            k=20,intensity=FALSE, sparse = TRUE)
#' plot(densityOut)
#' homerange(densityOut, percent = 0.95)
#' 
#' 
#' @import utils
#' @import graphics
#' @import stats
#' @importFrom spdep dnearneigh
#' @export
formLattice <-
function(nodeFillingOutput){
   if(class(nodeFillingOutput)!="nodeFillingOutput"){
       stop("Should be the output from the functions 
            nodeFilling or removeHole")
     }
   nodes <- nodeFillingOutput$nodes
   poly <- nodeFillingOutput$poly
   node_spacing <- nodeFillingOutput$node_spacing
   latt <- spdep::dnearneigh(nodes, node_spacing*0.5, node_spacing*1.5)
   formLatticeOutput <- list(EW_locs = nodeFillingOutput$EW_locs,
                              NS_locs = nodeFillingOutput$NS_locs,
                              nodes = nodes,
                              poly = poly,
                              latt = latt,
                              hole_list = nodeFillingOutput$hole_list)
   class(formLatticeOutput) <- "formLatticeOutput"
   return(formLatticeOutput)
}



