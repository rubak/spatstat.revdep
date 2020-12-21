#' UBC crossvalidation for the lattice-based density estimator.
#'
#' A function to perform crossvalidation to determine the smoothing
#' parameter for the lattice-based density estimator.  It minimizes the
#' UCV criterion.
#' 
#' @details The function computes the k-step diffusion \eqn{p_k = T^kp_0}, then computes the
#' Unbiased CrossValidation (UCV) criterion of Sain, Baggerly and Scott (1994).
#' This function can compute the UCV using either full matrix methods or sparse 
#' (default) matrix methods.  The latter are almost always much faster, though it 
#' is possible that if the number of points in the point pattern is large compared
#' to the number of nodes (an unlikely circumstance) that the full matrix method
#' would be quicker.  The sparse matrix approach typically uses less memory.  The
#' paper by Barry and McIntyre (2010) shows the approximation to the UCV used
#' in this approach.
#' @param formLatticeOutput An object from formLattice or editLattice.
#' @param PointPattern A matrix or data frame of locations.
#' @param M The maximum probability that the random walk will move.
#' @param max_steps The maximum number of steps attempted.
#' @param sparse Whether spare matrix computations used.
#' @return 
#' \itemize{
#'   \item ucv The value of the goodness-of-fit statistic.
#'   \item k The number of steps.
#' }
#' @references 
#' Crossvalidation of Multivariate Densities.
#' Stephan R. Sain, Keith A. Baggerly, David W. Scott; Journal of the American 
#' Statistical Association, Vol. 89 (1994) 807-817
#' 
#' @references Julie McIntyre, Ronald P. Barry (2018)  A Lattice-Based 
#' Smoother for Regions with Irregular Boundaries and Holes.  
#' Journal of Computational and Graphical Statistics.  In Press.
#' 
#' @examples 
#' plot.new()
#' data(polygon1)
#' #
#' nodeFillingOutput <- nodeFilling(poly=polygon1, node_spacing=0.02)
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
#' out <- crossvalDensity(formLatticeOutput,PointPattern=Pointdata, 
#'        M=0.5,max_steps = 70)
#' #
#' densityOut <- createDensity(formLatticeOutput,PointPattern=Pointdata, 
#'                           k=out$k,intensity=FALSE, sparse = TRUE)
#'plot(densityOut)
#' #
#' homerange(densityOut, percent = 0.95) 
#' @import utils
#' @import graphics
#' @import stats
#' @export
crossvalDensity <- 
function(formLatticeOutput, PointPattern, M=0.5, 
         max_steps = 200, sparse=TRUE)
  {
#
#
#
  if(class(formLatticeOutput)!="formLatticeOutput"){
       stop("Should be the output from the function formLattice")
    }
  if((M == 0)|(M == 1)){warning("Setting M to zero or one is ill-advised")
    }
  PointPattern <- as.matrix(PointPattern)
  addObsOutput <- addObservations(formLatticeOutput, PointPattern)
  p0 <- addObsOutput$init_prob
  which_nodes <- addObsOutput$which_nodes
  NN <- length(p0)
  n <- length(PointPattern[,1])
  poly_area <- areaRegion(formLatticeOutput)
#
if(sparse){
  #
  hold_del_prob <- matrix(nrow=NN, ncol=n, NA)
  for (i in 1:n){
   temp <- addObservations(formLatticeOutput, PointPattern[-i,,drop=FALSE])
   hold_del_prob[,i] <- temp$init_prob
  }
  T <- makeTmatrix(formLatticeOutput, M = M, sparse=TRUE)
  Tk_init <-  p0
  first_term <- rep(NA, max_steps)
  second_term <- rep(NA, max_steps)
  for(k in 1:max_steps){
    Tk_init <- T%*%Tk_init
    first_term[k] <- (NN/poly_area)*sum(Tk_init*Tk_init)
    hold_del_prob <- T%*%hold_del_prob
    temp <- sum(hold_del_prob[cbind(which_nodes, 1:n)])
    second_term[k] <- (NN/poly_area)*(2/n)*temp
} 
#
} else {
#
  T <- makeTmatrix(formLatticeOutput, M = 0.5, sparse=FALSE)
  Tnew <- diag(x=1, nrow=NN, ncol=NN)
  first_term <- rep(NA, max_steps)
  second_term <- rep(NA, max_steps)
  for(k in 1:max_steps){
    Tnew <- Tnew%*%T
    probs <- Tnew%*%p0
    first_term[k] <- (NN/poly_area)*sum(probs^2)
    second_term[k] <- (NN/poly_area)*(2/n)*((n^2/(n-1))*p0%*%probs - (n/(n-1))*sum(diag(Tnew)*p0))
  } 
 }
  ucv <- first_term - second_term
  k <- which.min(ucv)
  out <- list(ucv=ucv, k=k)
  return(out)
  }