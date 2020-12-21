#' Input data for Nonparametric Regression smoothing.
#' 
#' This function takes a \code{formLatticeOutput} object, which encodes a region
#' possibly with and irregular boundary and holes.  This and a matrix of
#' locations where a response variable has been measured, and a vector of
#' the responses, is used to create an initial distribution for use in the
#' non-parametric regression function \code{createNparReg}.  If, instead, you
#' have a point process and wish to produce a density estimate, you should use
#' the function \code{addObservations}.
#' 
#' @param formLatticeOutput An object from the functions formLattice or
#' editLattice.
#' @param Z A vector of response variable values.
#' @param locations A two-column matrix or data frame of data locations.
#' 
#' #' @references Ronald P. Barry, Julie McIntyre.  Estimating animal densities and home
#' range in regions with irregular boundaries and holes:  A lattice-based
#' alternative to the kernel density estimator.
#' Ecological Modelling 222 (2011)  1666-1672.
#' 
#' @references Julie McIntyre, Ronald P. Barry (2018)  A Lattice-Based 
#' Smoother for Regions with Irregular Boundaries and Holes.  
#' Journal of Computational and Graphical Statistics.  In Press.
#' 
#' @import utils
#' @import graphics
#' @import stats
#' @import spatstat
#' @import sp
#' @export
addQuantVar <-
function(formLatticeOutput, Z, locations){
  #
  #
  if(class(formLatticeOutput) != "formLatticeOutput"){
       stop("Should be the output from the function formLattice")
    }
  nodes <- formLatticeOutput$nodes
  n_observ <- nrow(locations)
  if(length(Z) != n_observ){stop("The number of rows in the argument
  locations should be the same as the length of Z")
    }
  n_nodes <- nrow(nodes)
  Z <-  as.vector(Z)
#
#  Here we have to create ppp objects to use spatstat functions
#  to find nearest latts in nodes from locations
#
  temp <- sp::bbox(rbind(locations,nodes))
  bound_vect <- c(temp[1,1], temp[1,2], temp[2,1], temp[2,2])
  X <- spatstat.geom::as.ppp(locations, W=bound_vect)
  Y <- spatstat.geom::as.ppp(nodes, W=bound_vect)
  closest <- spatstat.geom::nncross(X,Y)$which
#
#  The output will be a vector that gives an initial prob
#  at each node, depending on number of corresponding
#  locations
#
  temp <- addObservations(formLatticeOutput, observations=locations)
  sums <- function(x){sum(Z[closest==x])}
  out <- list(init_quantvar = as.numeric(lapply(1:n_nodes, FUN=sums))/n_observ,
             init_prob = temp$init_prob,
             which_nodes = closest)
  class(out) <- "addQuantVarOutput"
  return(out)
}
