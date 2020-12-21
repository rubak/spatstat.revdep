#' Plot the non-parametric regression surface.
#' 
#' Takes as input a NparRegOut object from the 
#' function createNparReg. This plotting function 
#' makes a countour plot of the non-parametric 
#' regression prediction surface.
#' 
#' @references Julie McIntyre, Ronald P. Barry (2018)  A Lattice-Based 
#' Smoother for Regions with Irregular Boundaries and Holes.  
#' Journal of Computational and Graphical Statistics.  In Press.
#' @param x An object of type NparRegOut returned by createNparReg.
#' @param \dots Other arguments to be passed to functions plot, 
#' points, lines.
#' @author Ronald P. Barry
#' @examples 
#' data(nparExample)
#' attach(nparExample)
#' plot.new()
#' #  Simulate a response variable
#' index1 <- (grid2[,2]<0.8)|(grid2[,1]>0.6)
#' Z <- rep(NA,length(grid2[,1]))
#' n1 <- sum(index1)
#' n2 <- sum(!index1)
#' Z[index1] <- 3*grid2[index1,1] + 4 + rnorm(n1,0,sd=0.4)
#' Z[!index1] <- -2*grid2[!index1,1] + 4 + rnorm(n2,0,sd=0.4)
#' #
#' plot(polygon2,type="n")
#' polygon(polygon2)
#' points(grid2,pch=19,cex=0.5,xlim=c(-0.1,1))
#' text(grid2,labels=round(Z,1),pos=4,cex=0.5)
#' #  Following is the generation of the nonparametric
#' #  regression prediction surface.
#' nodeFillingOutput <- nodeFilling(poly=polygon2, node_spacing=0.025)
#' plot(nodeFillingOutput)
#' formLatticeOutput <- formLattice(nodeFillingOutput)
#' plot(formLatticeOutput)
#' NparRegOut <- createNparReg(formLatticeOutput,Z,PointPattern=grid2,k=2)
#' plot(NparRegOut)
#' @import utils
#' @import graphics
#' @import stats
#' @export
#' @method plot NparRegOut
plot.NparRegOut <-
function(x,...){
#
  if(class(x) != "NparRegOut"){
    stop("This function only plots objects of type NparRegOut")
    }
  EW_locs <- as.vector(x$EW_locs)
  NS_locs <- as.vector(x$NS_locs)
  NparRegMap <- as.vector(x$NparRegMap)
  boundaryPoly <- as.matrix(x$boundaryPoly)
  hole_list <- x$hole_list
  big_matrix <- matrix(nrow=length(EW_locs),ncol=length(NS_locs),NparRegMap)
  contour(x=EW_locs,y=NS_locs,z=big_matrix,...)
  lines(rbind(boundaryPoly,boundaryPoly[1,]),...)
  n_hole <- length(hole_list)
  for(i in 1:n_hole){
    lines(rbind(hole_list[[i]],hole_list[[i]][1,]),...)
    }
  }
  
