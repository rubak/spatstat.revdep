#' Performs nonparametric regression on irregular regions.
#' 
#' This function takes the lattice from formLattice (which fills the 
#' region of interest) along with the list of responses and their locations,
#' and creates a prediction surface.  The approach is kernel non-parametric
#' regression with the kernels created by a k-step diffusion on the
#' lattice about each location where a response was collected.
#' 
#' We denote by \eqn{K_{ik}(s)} the kernel obtained by assigning the node
#' nearest to the ith response and then running a k-step diffusion on 
#' the lattice and evaluating the resulting density at location s.  
#' Then the estimator \eqn{\hat{f}(s) = (\sum_i K_{ik}(s)*Z_i)/\sum_i K_{ik}(s)}
#' which is the traditional kernal regression estimator with diffusion 
#' kernels.  This approach leads to a non-parametric regression that
#' respects the boundaries of the polygonal region.  The construction of the
#' kernels is detailed in Barry and McIntyre (2011).  Using kernels to 
#' perform nonparametric regression is described in many publications, including
#' Wasserman (2006).
#' @section Variance Estimation:  We use the variance estimator \eqn{\sum e_{i,-i}^2/n}, where \eqn{e_{i,-i}} is the ith deleted residual.
#' 
#' @references
#' Larry Wasserman.  All of Nonparametric Statistics.  Springer Science + 
#' Business Media, Inc. N.Y. 2006.
#' 
#' #' @references Julie McIntyre, Ronald P. Barry (2018)  A Lattice-Based 
#' Smoother for Regions with Irregular Boundaries and Holes.  
#' Journal of Computational and Graphical Statistics.  In Press.
#' 
#' @param formLatticeOutput An object returned by formLattice or editLattice.
#' @param Z Vector of responses to be smoothed.
#' @param PointPattern A 2 column matrix or data frame of locations.
#' @param M The maximum probability that the random walk will move.
#' @param k Number of steps.
#'   
#' @examples
#' data(nparExample)
#' attach(nparExample)
#' plot.new()
#' #  Simulate a response variable
#' index1 = (grid2[,2]<0.8)|(grid2[,1]>0.6)
#' Z = rep(NA,length(grid2[,1]))
#' n1 = sum(index1)
#' n2 = sum(!index1)
#' Z[index1] <- 3*grid2[index1,1] + 4 + rnorm(n1,0,sd=0.4)
#' Z[!index1] <- -2*grid2[!index1,1] + 4 + rnorm(n2,0,sd=0.4)
#' #
#' coords=rbind(polygon2,polygon2[1,])
#' plot(coords,type="l")
#' points(grid2,pch=19,cex=0.5,xlim=c(-0.1,1))
#' text(grid2,labels=round(Z,1),pos=4,cex=0.5)
#' #
#' nodeFillingOutput <- nodeFilling(poly=polygon2, node_spacing=0.025)
#' plot(nodeFillingOutput)
#' formLatticeOutput <- formLattice(nodeFillingOutput)
#' NparRegOut <- createNparReg(formLatticeOutput,Z,PointPattern=grid2,k=2)
#' plot(NparRegOut)
#'
#' @import utils
#' @import graphics
#' @import stats
#' @export
createNparReg <-
function(formLatticeOutput, Z, 
         PointPattern=NULL, M=0.5, k)
{
#
#  This function only uses sparse matrices
#  It takes as input the formLatticeOutput and a column of responses
#  called Z.  Finally, k is the smoothing parameter and is the
#  number of steps in the random walk.
#
  if(class(formLatticeOutput) != "formLatticeOutput"){
       stop("Should be the output from the function formLattice")
    }
  if((M == 0)|(M == 1)){warning("Setting M to zero or one is ill-advised")
    }
  if(!is.null(PointPattern)){PointPattern <- as.matrix(PointPattern)}
  addQuantVarOutput <- addQuantVar(formLatticeOutput = formLatticeOutput, 
                       Z=Z,
                       locations = PointPattern)
  p0 <- addQuantVarOutput$init_prob
  Z0 <- addQuantVarOutput$init_quantvar
  which_nodes <- addQuantVarOutput$which_nodes
  T<-makeTmatrix(formLatticeOutput=formLatticeOutput, M=M, sparse=TRUE)
  EW_locs <- formLatticeOutput$EW_locs
  NS_locs <- formLatticeOutput$NS_locs
  nodes <- formLatticeOutput$nodes
  pk <- Tkp(T, k=k, p=p0)
  zk <- Tkp(T, k=k, p=Z0)
  N <- length(NS_locs)*length(EW_locs)
  long <- rep(NA,N)
  long[as.numeric(rownames(nodes))] <- zk/pk
  long[is.nan(long)] <- mean(Z)
  #  Compute an estimate of variability
  deleted_residuals <- 
    deletedResid(formLatticeOutput,Z,PointPattern,M=M,k=k)$deleted_residuals
  sigma2 <- mean(deleted_residuals^2)
  NparRegOut <- list(EW_locs = formLatticeOutput$EW_locs,
                    NS_locs = formLatticeOutput$NS_locs,
                    nodes = formLatticeOutput$nodes,
                    boundaryPoly = formLatticeOutput$poly,
                    hole_list = formLatticeOutput$hole_list,
                    PointPattern = PointPattern,
                    which_nodes = which_nodes,
                    NparRegNum = zk,
                    NparRegDenom = pk,
                    NparRegMap = long,
                    sigma2 = sigma2)
  class(NparRegOut) <- "NparRegOut"
  return(NparRegOut)
}


