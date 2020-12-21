#' Deleted residuls for non-parametric regression.
#' 
#' Computes deleted residuals for the lattice-based 
#' non-parametric regression estimator.
#' 
#' @param formLatticeOutput An object from formLattice or editLattice.
#' @param Z Vector of response values.
#' @param PointPattern 2-column matrix or data frame of locations.
#' @param M Maximum probability that the random walk moves.
#' @param k Number of steps in random walk.
#' @return A vector of deleted residuals.
#' @author Ronald P. Barry
#' @references 
#' Ronald P. Barry, Julie McIntyre. Estimating 
#' animal densities and home range in regions 
#' with irregular boundaries and holes: A 
#' lattice-based alternative to the kernel 
#' density estimator. Ecological Modelling 222 (2011) 1666-1672.
#' 
#' @references Julie McIntyre, Ronald P. Barry (2018)  A Lattice-Based 
#' Smoother for Regions with Irregular Boundaries and Holes.  
#' Journal of Computational and Graphical Statistics.  In Press.
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
#' #
#' nodeFillingOutput <- nodeFilling(poly=polygon2, node_spacing=0.025)
#' plot(nodeFillingOutput)
#' formLatticeOutput <- formLattice(nodeFillingOutput)
#' plot(formLatticeOutput)
#' hold <- crossvalNparReg(formLatticeOutput,Z,
#'                        PointPattern=grid2,M=0.5,max_steps = 40)
#' deletedResid(formLatticeOutput,Z,
#'              PointPattern=grid2,M=0.5,k=hold$k)
#' 
#' @import utils
#' @import graphics
#' @import stats
#' @export
deletedResid <- 
function(formLatticeOutput,Z,PointPattern,M=0.5,k){
#
#                                                            
#
  if(class(formLatticeOutput)!="formLatticeOutput"){
       stop("Should be the output from the function formLattice")}
  if((M==0)|(M==1)){
               warning("Setting M to zero or one is ill-advised")}
  PointPattern <- as.matrix(PointPattern)
  addQuantVarOut <- addQuantVar(formLatticeOutput, Z=Z, 
      locations = PointPattern)
  #  Despite being called pk and Zk, these are initially the
  #  starting values.  Later they will be updated k times.
  pk <- addQuantVarOut$init_prob
  Zk <- addQuantVarOut$init_quantvar
  which_nodes <- addQuantVarOut$which_nodes
  NN <- length(pk)
  n <- length(PointPattern[,1])
  del_prob <- matrix(nrow=NN, ncol=n, NA)
  del_Z <- matrix(nrow=NN, ncol=n, NA)
  #  Create initial probability and quantitative variable
  #  vectors for each deleted observation.
  for (i in 1:n){
    temp <- addQuantVar(formLatticeOutput, Z=Z[-i], 
        locations = PointPattern[-i,,drop=FALSE])
    del_prob[,i] <- temp$init_prob
    del_Z[,i] <- temp$init_quantvar
  }
  #  Construct the overall transition matrix, then compute
  #  k steps in the random walk for original data and all
  #  deleted cases.
  T <- makeTmatrix(formLatticeOutput, M = M, sparse=TRUE)
  for(k in 1:k){
    pk <- T%*%pk
    zk <- T%*%Zk
    del_prob <- T%*%del_prob
    del_Z <- T%*%del_Z
  }
#   fitted values are computed at locations with observations
    fitted <- (zk/pk)[which_nodes,]
    fitted[is.nan(fitted)] <- mean(Z)
    residuals <- Z - fitted
    #  Extract the deleted predictions from the matrix of predictions
    #  under all deleted models.  Note that some will be NAN, due to
    #  division by zero.  These will be replaced by yi minus the mean
    #  of all other observations.
    del_predictions <- diag((del_Z/del_prob)[which_nodes,])
    nnan <- sum(is.nan(del_predictions))
    del_predictions[is.nan(del_predictions)] <- (rep(sum(Z),nnan) - 
             Z[is.nan(del_predictions)])/(n-1)
    deleted_residuals <- Z - del_predictions
   return(list(deleted_residuals=deleted_residuals,residuals=residuals)) 
}











