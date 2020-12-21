#' Input observations for use in the lattice-based density estimator
#'
#' This function takes a \code{formLatticeOutput} object, which encodes a region
#' possibly with and irregular boundary and holes, and adds point process
#' observations.  The observations should be in the form of a matrix or
#' data frame.  \code{addObservations} is used when the aim is to produce
#' a density map from a point process.  If, instead, you wish to produce
#' a nonparametric regression surface given responses and their locations, you
#' should use \code{addQuantVar} instead.
#'  
#' @details
#' Every node in the \code{formLatticOutput} object is assigned an initial density
#' that is equal to the fraction of all observations that are nearest to that
#' node.  Note that this means observations can be outside the boundary of the
#' region of interest - they will just be associated with the nearest node inside
#' the region.  The function returns a vector equal in length to the number of
#' nodes that has the initial density for each node.  This vector corresponds to
#' \eqn{p_0}, the inital probability vector as in Barry and McIntyre (2011).
#
#' @param formLatticeOutput  An object returned by formLattice or
#' editLattice.
#' @param observations A matrix or data frame with two columns.
#'
#' @return a list with two elements.
#' \itemize{
#' \item{init_prob}{Numerical vector with the initial probability distribution}
#' \item{which_nodes}{vector of nodes to which observations were assigned}
#' }
#' @references Ronald P. Barry, Julie McIntyre.  Estimating animal densities and home
#' range in regions with irregular boundaries and holes:  A lattice-based
#' alternative to the kernel density estimator.
#' Ecological Modelling 222 (2011)  1666-1672.
#'
#' @author Ronald P. Barry
#
#' @examples
#' plot.new()
#' data(polygon1)
#' #
#' nodeFillingOutput <- nodeFilling(poly=polygon1, node_spacing=0.01)
#' plot(nodeFillingOutput)
#' formLatticeOutput <- formLattice(nodeFillingOutput)
#' plot(formLatticeOutput)
#' #
#' Pointdata <- splancs::csr(polygon1,30)
#' Pointdata <- Pointdata[Pointdata[,1]<0.5,]
#' colnames(Pointdata) <- c("x","y")
#' plot(polygon1,type="n")
#' polygon(polygon1)
#' points(Pointdata,pch=19)
#' #
#' densityOut <- createDensity(formLatticeOutput,PointPattern=Pointdata,
#'                            k=40,intensity=FALSE, sparse = TRUE)
#' plot(densityOut)
#' @import utils
#' @import graphics
#' @import stats
#' @import spatstat
#' @import sp
#' @export
addObservations <-
function(formLatticeOutput, observations){
  #
  if(class(formLatticeOutput) != "formLatticeOutput"){
       stop("Should be the output from the function formLattice")
    }
  nodes <- formLatticeOutput$nodes
  n_observ <- nrow(observations)
  n_nodes <- nrow(nodes)
#
#  Here we have to create ppp objects to use spatstat functions
#  to find nearest latts in nodes from observations
#
  temp <- sp::bbox(rbind(observations, nodes))
  bound_vect <- c(temp[1,1],temp[1,2],temp[2,1],temp[2,2])
  X <- spatstat::as.ppp(observations, W=bound_vect)
  Y <- spatstat::as.ppp(nodes, W=bound_vect)
  closest <- spatstat::nncross(X,Y)$which
#
#  The output will be a vector that gives an initial prob
#  at each node, depending on number of corresponding
#  observations
#
  out <- list(init_prob = tabulate(closest,nbins=n_nodes)/n_observ,
             which_nodes = closest)
  class(out) <- "initProbObject"
  return(out)
}

