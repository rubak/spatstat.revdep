#' Spatial variance for the regression smoother.
#' 
#' Computes the variance at each location for the 
#' non-parametric regression estimator.
#' 
#' @details  \code{varianceMap} computes an estimated variance at each
#' node in the lattice, output in a form for mapping with contour.
#' The approach is the Nadaraya-Watson kernel variance estimator:
#' \eqn{s^2\sum K^2(si,s0)/(\sum K(si,s0))^2}.  It's important to note
#' that this should not be overused as a prediction error, as kernel
#' estimators are not unbiased.
#' 
#' @param formLatticeOutput An object from formLattice or editLattice.
#' @param Z Vector of response values.
#' @param PointPattern 2-column matrix or data frame of locations.
#' @param M Maximum probability that the random walk moves.
#' @param k Number of steps in random walk.
#' @return VarianceMapOut object
#' \itemize{
#'   \item EW_locs EW coordinates for use by contour
#'   \item NS_locs NS coordinates for use by contour
#'   \item boundaryPoly vertices of the boundary
#'   \item hole_list list of polygonal hole boundaries, if any.
#'   \item SE_map_grid estimated standard error at each location}
#'   
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
#'
#' polygon(polygon2)
#' points(grid2,pch=19,cex=0.5,xlim=c(-0.1,1))
#' text(grid2,labels=round(Z,1),pos=4,cex=0.5)
#' #
#' nodeFillingOutput <- nodeFilling(poly=polygon2, node_spacing=0.025)
#' plot(nodeFillingOutput)
#' formLatticeOutput <- formLattice(nodeFillingOutput)
#' plot(formLatticeOutput)
#' hold <- crossvalNparReg(formLatticeOutput,Z,
#'                        PointPattern=grid2,M=0.5,max_steps = 20)
#' var_map <- varianceMap(formLatticeOutput,Z,
#'              PointPattern=grid2,M=0.5,k=10)
#' plot(var_map)
#' 
#' @import utils
#' @import graphics
#' @import stats
#' @export
varianceMap <- 
  function(formLatticeOutput,Z,PointPattern,M=0.5,k){
    #
    #                                                            
    #
    if(class(formLatticeOutput)!="formLatticeOutput"){
      stop("Should be the output from the function formLattice")
      }
    if((M==0)|(M==1)){
      warning("Setting M to zero or one is ill-advised")}
    PointPattern <- as.matrix(PointPattern)
    NS_locs <- formLatticeOutput$NS_locs
    EW_locs <- formLatticeOutput$EW_locs
    nodes <- formLatticeOutput$nodes
    boundaryPoly <- formLatticeOutput$poly
    hole_list <- formLatticeOutput$hole_list
    #
    #  Get the crossvalidated sum of squared error
    #
    crossval_sum_square <- 
                    crossvalNparReg(formLatticeOutput=formLatticeOutput, Z=Z, 
                    PointPattern=PointPattern,M=M,max_steps = k)$SumSq[k]
    #
    NN <- length(nodes[,1])
    n <- length(PointPattern[,1])
    var_prob <- matrix(nrow=NN, ncol=n, 0)
    #  Create initial probability and quantitative variable
    #  vectors for each observation.
    for (i in 1:n){
      temp <- addQuantVar(formLatticeOutput, Z=Z[i], 
                          locations = PointPattern[i,, drop=FALSE])
      var_prob[,i] <- temp$init_prob
    }
    #  Construct the overall transition matrix, then compute
    #  k steps in the random walk for original data and all
    #  deleted cases.
    T <- makeTmatrix(formLatticeOutput, M = M, sparse=TRUE)
    #
    for(k in 1:k){
      var_prob <- T%*%var_prob
    }
    #
    #  var_prob contains the weights from n kernels, K(si,s0)
    #  crossval_sum_square has the crossvalidated SSE
    sigma2_hat <- crossval_sum_square/n
    numerator <- rowSums(var_prob^2)*sigma2_hat
    no_weights <- numerator <= 0
    denominator <- rowSums(var_prob)^2
    SE_est <- sqrt(numerator/denominator)
    SE_est[no_weights] <- NA
    SE_map_size <- length(NS_locs)*length(EW_locs)
    SE_map <- rep(NA, SE_map_size)
    SE_map[as.numeric(rownames(nodes))] <- SE_est
  varianceMapOut = list(EW_locs = EW_locs,
    NS_locs = NS_locs,
    boundaryPoly = boundaryPoly,
    hole_list = hole_list,
    SE_map = SE_map
    )
class(varianceMapOut) = "varianceMapOut"
return(varianceMapOut)
  }
