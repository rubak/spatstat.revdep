#' Plot the lattice.
#' 
#' This function plots the boundary, holes, 
#' nodes and neighbor lattice for the lattice 
#' based density or regression estimators. The 
#' plot can be examined to determine whether 
#' the lattice of connected nodes fills the region. 
#' If some nodes are connected when they should 
#' not be, or are disconnected when they should be 
#' connected, use editLattice to add or remove neighbor links.
#' 
#' @param x An object of type formLatticeOutput returned 
#' by either formLattice or editLattice.
#' @param \dots Other arguments to be passed to functions 
#' plot, points, lines.
#' @author Ronald P. Barry
#' @examples 
#' plot.new()
#' data(polygon1)
#' nodeFillingOutput <- nodeFilling(poly=polygon1, node_spacing=0.015)
#' plot(nodeFillingOutput)
#' formLatticeOutput <- formLattice(nodeFillingOutput)
#' plot(formLatticeOutput)
#' 
#' @import utils
#' @import graphics
#' @import stats
#' @importFrom spdep card
#' @export
#' @method plot formLatticeOutput
plot.formLatticeOutput <-
function(x,...){
  #  
  #  This function plots the lattice.  If there are
  #  unwanted lattice links, you may edit the lattice using
  #  editFromLatticeOutput
  #
  if(class(x)!="formLatticeOutput"){
       stop("Should be the output from the function formLattice")}
  nodes <- x$nodes
  poly <- x$poly
  latt <- x$latt
  xp_min <- min(poly[,1])
  xp_max <- max(poly[,1])
  yp_min <- min(poly[,2])
  yp_max <- max(poly[,2])
  old_plt <- par()$plt
  on.exit(par(plt=old_plt))
  par(plt=c(0,1,0,1))
  plot(coords=nodes,latt,xlim=c(xp_min,xp_max), ylim= c(yp_min, yp_max),
    cex=0.1,...)
  lines(rbind(poly,poly[1,]),lwd=1,...)
  if(min(spdep::card(latt))==0){
      points(nodes[spdep::card(latt)==0,],pch=19,...)
  }
}

