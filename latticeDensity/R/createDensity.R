#' Generates a density using random walks on a lattice.
#' 
#' Given a lattice and a point pattern of observations, 
#' createDensity starts random walks at each observation. 
#' k steps are taken and the output is a densityOut object, 
#' which can be used to plot a density estimate. If you wish 
#' to perform non-parametric regression, you should use the 
#' functions addQuantVar and createNparReg instead.
#' 
#' @details
#' We start with an initial probability density p0 where the 
#' ith entry in p0 is the fraction of the point pattern that 
#' is nearest to the ith node. This is the empirical density 
#' function with no smoothing. If T is the transition matrix, 
#' and given a number of steps in the diffusion, T k p0 is the 
#' probability density of the diffusion after k steps. This is 
#' the major output of this function, along with information 
#' needed to produce a plot, including the polygons for the 
#' boundary and holes, and a vector of NS coordinates and EW 
#' coordinates used by the contour function. All of the necessary 
#' information for plotting is bundled in the object of class 
#' densityOutLBDE. Details of this process can be found in Barry 
#' and McIntyre (2011).
#' @param formLatticeOutput An object from formLattice or editLattice.
#' @param PointPattern A 2-column matrix or data frame of locations.
#' @param M Maximum probability of random walk moving.
#' @param k The smoothing parameter (number of steps).
#' @param intensity Plot an intensity vs a density.
#' @param sparse If TRUE, matrix computations are sparse.
#' @return An object of type densityOut
#' \itemize{
#'   \item EW_locs A vector of EW coordinates of nodes.
#'   \item NS_locs A vector of NS coordinates of nodes.
#'   \item boundaryPoly The boundary of the region (two-columns).
#'   \item hole_list A list of polygonal holes in the region.
#'   \item PointPattern A 2-column matrix of observations.
#'   \item probs The probability distribution over the nodes.
#'   \item densityLBDE Density in a form for making a contour map.
#'   \item area The area of the region, with holes removed.}
#' @author Ronald P. Barry
#' @references
#' Ronald P. Barry, Julie McIntyre. Estimating animal 
#' densities and home range in regions with irregular 
#' boundaries and holes: A lattice-based alternative to 
#' the kernel density estimator. Ecological Modelling 
#' 222 (2011) 1666-1672.
#' @examples 
#' plot.new()
#' data(polygon1)
#' nodeFillingOutput <- nodeFilling(poly=polygon1, node_spacing=0.02)
#' plot(nodeFillingOutput)
#' formLatticeOutput <- formLattice(nodeFillingOutput)
#' plot(formLatticeOutput)
#' Pointdata <- splancs::csr(polygon1,75)
#' Pointdata <- Pointdata[Pointdata[,1]<0.5,]
#' plot(polygon1,type="n")
#' polygon(polygon1)
#' points(Pointdata,pch=19)
#' out <- crossvalDensity(formLatticeOutput,PointPattern=Pointdata,
#'                       M=0.5,max_steps = 35)
#' densityOut <- createDensity(formLatticeOutput,
#'                            PointPattern=Pointdata, 
#'                            k=out$k,intensity=FALSE, sparse = TRUE)
#' plot(densityOut)
#' homerange(densityOut, percent = 0.95)
#' 
#' 
#' @import utils
#' @import graphics
#' @import stats
#' @export
createDensity <-
function(formLatticeOutput, PointPattern=NULL, M=0.5, k,
         intensity=FALSE, 
         sparse=TRUE)
{
#
  if(class(formLatticeOutput)!="formLatticeOutput"){
       stop("Should be the output from the function formLattice")
    }
  if((M == 0)|(M == 1)){warning("Setting M to zero or one is ill-advised")
    }
  init_prob <- addObservations(formLatticeOutput, PointPattern)
  p0 <- as.vector(init_prob$init_prob)
  if(!is.null(PointPattern)){PointPattern <- as.matrix(PointPattern)}
  poly_area <- areaRegion(formLatticeOutput)
  n <- length(PointPattern[,1])
  T<-makeTmatrix(formLatticeOutput, M=M, sparse=sparse)
  EW_locs <- formLatticeOutput$EW_locs
  NS_locs <- formLatticeOutput$NS_locs
  nodes <- formLatticeOutput$nodes
  z <- Tkp(T, k=k, p=p0)
  N <- length(NS_locs)*length(EW_locs)
  long <- rep(NA,N)
  if(intensity){
  long[as.numeric(rownames(nodes))] <-
       n*z*length(z)/poly_area
  }else{
  long[as.numeric(rownames(nodes))] <-
       z*length(z)/poly_area
  }
  densityOut     <- list(EW_locs = formLatticeOutput$EW_locs,
                    NS_locs = formLatticeOutput$NS_locs,
                    nodes = formLatticeOutput$nodes,
                    boundaryPoly = formLatticeOutput$poly,
                    hole_list = formLatticeOutput$hole_list,
                    PointPattern = PointPattern,
                    probs = z,
                    densityOut = long,
                    area = poly_area)
  class(densityOut) <- "densityOut"
  return(densityOut)
}


