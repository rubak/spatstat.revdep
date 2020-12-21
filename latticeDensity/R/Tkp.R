#' Compute the vector T^k*p
#' 
#' T is the transition matrix of the 
#' random walk on the lattice. By multiplying 
#' by the probability density p at time t, 
#' you get the probability density at time t+1. 
#' Thus, to get the probability density after k 
#' steps, pk, compute pk = Tkp1. This application 
#' of finite Markov processes is described in Barry 
#' and McIntyre (2011).
#' @param T Transition matrix returned by makeTmatrix.
#' @param k The number of steps in the diffusion.
#' @param p A numerical vector of length equal to the 
#' number of nodes, of initial probabilities.
#' @author Ronald P. Barry
#' @references Ronald P. Barry, Julie McIntyre. Estimating 
#' animal densities and home range in regions with irregular 
#' boundaries and holes: A lattice-based alternative to the 
#' kernel density estimator. Ecological Modelling 222 (2011) 
#' 1666-1672.
#' @examples 
#' plot.new()
#' data(polygon1)
#' nodeFillingOutput <- nodeFilling(poly=polygon1, node_spacing=0.015)
#' formLatticeOutput <- formLattice(nodeFillingOutput)
#' Pointdata <- splancs::csr(polygon1,75)
#' Pointdata <- Pointdata[Pointdata[,1]<0.5, ]
#' init_prob <- addObservations(formLatticeOutput, Pointdata)
#' T <- makeTmatrix(formLatticeOutput, M = 0.5, sparse=TRUE)
#' p10 <- Tkp(T, k=10, p=init_prob$init_prob)
#' head(cbind(init_prob$init_prob, p10))
#' 
#' @import utils
#' @import graphics
#' @import stats
#' @importFrom spam is.spam
#' @importFrom spam diag
#' @export
Tkp <-
function(T,k,p){
#  p should be output from 
  k <- as.numeric(k)
  if (class(p)=="initProbObject"){
    p <- as.vector(p$init_prob)
    } else {
    p <- as.vector(p)
    }
if(spam::is.spam(T)){
  # sparse matrix computations
  TkpOut <- p
  for (i in 1:k){
    TkpOut <- T%*%TkpOut
  }
} else {
  # full matrix computations
  T <- as.matrix(T)
  TkpOut <- p
  for (i in 1:k){
    TkpOut <- T%*%TkpOut
  }
}
return(TkpOut)
}

