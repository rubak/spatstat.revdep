#' Plot the density map
#' 
#' Plots the boundary, all holes and the locations of 
#' all nodes along with the density contour map.
#' 
#' @param x An object of type densityOut returned by createDensity.
#' @param \dots Graphical parameters for the function contour.default.
#' 
#' @author Ronald P. Barry
#' @references Ronald P. Barry, Julie McIntyre. Estimating 
#' animal densities and home range in regions with irregular 
#' boundaries and holes: A lattice-based alternative to the 
#' kernel density estimator. Ecological Modelling 222 (2011) 
#' 1666-1672
#' @examples 
#' plot.new()
#' data(polygon1)
#' #
#' nodeFillingOutput <- nodeFilling(poly=polygon1, node_spacing=0.025)
#' plot(nodeFillingOutput)
#' formLatticeOutput <- formLattice(nodeFillingOutput)
#' plot(formLatticeOutput)
#' #
#' Pointdata <- splancs::csr(polygon1,75)
#' Pointdata <- Pointdata[Pointdata[,1]<0.5,]
#' plot(polygon1,type="n")
#' polygon(polygon1)
#' points(Pointdata,pch=19)
#' #
#' densityOut <- createDensity(formLatticeOutput,PointPattern=Pointdata,
#'                            k=55,intensity=FALSE, sparse = TRUE)
#' plot(densityOut)
#' #
#' homerange(densityOut, percent = 0.95)
#' 
#' @import utils
#' @import graphics
#' @import stats
#' @export
#' @method plot densityOut
plot.densityOut <-
function(x,...){
#
  if(class(x) != "densityOut"){
    stop("This function only plots objects of type densityOut")
    }
  EW_locs <- as.vector(x$EW_locs)
  NS_locs <- as.vector(x$NS_locs)
  densityOut <- as.vector(x$densityOut)
  boundaryPoly <- as.matrix(x$boundaryPoly)
  hole_list <- x$hole_list
  big_matrix <- matrix(nrow=length(EW_locs),ncol=length(NS_locs),densityOut)
  contour(x=EW_locs,y=NS_locs,z=big_matrix,...)
  lines(rbind(boundaryPoly,boundaryPoly[1,]),...)
  n_hole <- length(hole_list)
  for(i in 1:n_hole){
    lines(rbind(hole_list[[i]],hole_list[[i]][1,]),...)
    }
  }
  
