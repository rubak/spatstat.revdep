#' Crossvalidation for non-parametric regression.
#' 
#' Performs least-squares crossvalidation for the 
#' lattice-based non-parametric regression estimator.
#' 
#' @details 
#' For a given k, deleted residuals are computed for each of the observations.
#' The crossvalidation is based on minimization of the squares of the
#' deleted residuals.
#' 
#' @param formLatticeOutput An object from formLattice or editLattice.
#' @param Z Vector of response values to be smoothed.
#' @param PointPattern A 2-column matrix or data frame of locations.
#' @param M Maximum probability that the random walk moves.
#' @param max_steps Maximum number of steps attempted.
#' @return A list consisting of 
#' \itemize{
#'   \item SumSq Vector of crossvalidated sums of squares
#'   \item Number of steps that minimizes the crossvalidated SS.
#'   }
#' @author Ronald P. Barry
#' @references 
#' Ronald P. Barry, Julie McIntyre. Estimating animal 
#' densities and home range in regions with irregular 
#' boundaries and holes: A lattice-based alternative to 
#' the kernel density estimator. Ecological Modelling 
#' 222 (2011) 1666-1672. 
#' 
#'@references Julie McIntyre, Ronald P. Barry (2018)  A Lattice-Based 
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
#' NparRegOut <- createNparReg(formLatticeOutput,Z,PointPattern=grid2,k=hold$k)
#' plot(NparRegOut)
#' @import utils
#' @import graphics
#' @import stats
#' @export
crossvalNparReg <-
function(formLatticeOutput,Z,PointPattern,M=0.5,max_steps = 200){
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
  p0 <- addQuantVarOut$init_prob
  Z0 <- addQuantVarOut$init_quantvar
  which_nodes <- addQuantVarOut$which_nodes
  NN <- length(p0)
  n <- length(PointPattern[,1])
  del_prob <- matrix(nrow=NN, ncol=n, NA)
  del_Z <- matrix(nrow=NN, ncol=n, NA)
  for (i in 1:n){
    temp <- addQuantVar(formLatticeOutput, Z=Z[-i], 
        locations = PointPattern[-i,,drop=FALSE])
    del_prob[,i] <- temp$init_prob
    del_Z[,i] <- temp$init_quantvar
  }
  T <- makeTmatrix(formLatticeOutput,M = M, sparse=TRUE)
# # Tkinit <-  p0
# # TZinit <-  Z0
  SumSq <- rep(NA,max_steps)
  for(k in 1:max_steps){
  #  Tkinit <- T%*%Tkinit
  #  TZinit <- T%*%TZinit
    del_prob <- T%*%del_prob
    del_Z <- T%*%del_Z
    predictions <- diag((del_Z/del_prob)[which_nodes,])
    nnan <- sum(is.nan(predictions))
    predictions[is.nan(predictions)] <- (rep(sum(Z),nnan) - 
             Z[is.nan(predictions)])/(n-1)
    deleted_residuals <- predictions - Z
    SumSq[k] <- sum((deleted_residuals)^2)
   }
   plot(SumSq,type="l")
   k <- which.min(SumSq)
   out <- list(SumSq=SumSq,k=k)
   return(out) 
}











