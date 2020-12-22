#' Performs a sensitivity analysis on the parameters \code{k} and \code{n} that are provided to \code{DRHotspots_k_n}
#' 
#' Given a set of \code{ks} and \code{ns} parameters, this function allows the user to perform a sensitivity analysis on the parameters \code{k} and \code{n} by calling \code{DRHotspots_k_n} for each combination of \code{k} and \code{n}
#' 
#' @param X - A \code{lpp} object representing a marked point pattern lying on a road network (\code{linnet} object)
#' @param rel_probs - An object containing the relative probabilities of a specific type of event along the linear network contained in \code{X}, generated through the function \code{RelativeProbabilityNetwork}
#' @param ks - A numeric vector of possible values for the \code{k} parameter that is provided to \code{DRHotspots_k_n} 
#' @param ns - A numeric vector of possible values for the \code{n} parameter that is provided to \code{DRHotspots_k_n} 
#' @return A matrix providing the type-specific prediction accuracy index that corresponds to the set differential risk hotspots obtained for each value of \code{k} or \code{n} provided in \code{ks} and \code{ns}, respectively. A \code{NA} value in this matrix indicates that no differential risk hotspots are found for the corresponding combination of \code{k} and \code{n}
#' @examples 
#' library(DRHotNet)
#' library(spatstat.core)
#' library(spdep)
#' library(raster)
#' library(maptools)
#' \donttest{
#' rel_probs_rear_end <- RelativeProbabilityNetwork(X = SampleMarkedPattern,
#' lixel_length = 50,sigma = 100,mark = "Collision",category_mark = "Rear-end")
#' sensitivity_analysis <- Sensitivity_k_n(X = SampleMarkedPattern, rel_probs = rel_probs_rear_end, 
#' ks = c(1,2), ns = c(30,40))
#' }
#' @references Briz-Redon, A., Martinez-Ruiz, F., & Montes, F. (2019). Identification of differential risk hotspots for collision and vehicle type in a directed linear network. Accident Analysis & Prevention, 132, 105278.
#' @export
Sensitivity_k_n <- function(X,rel_probs,ks,ns){
  
  network=X$domain
  
  if (rel_probs$lixel_length!=F){
    network=lixellate(network,eps=rel_probs$lixel_length)
    #& project into the new network
    X_aux=lpp(cbind(X$data$x,X$data$y),network)
    marks(X_aux)=marks(X)
    X=X_aux
  }
  network_lix=X$domain
  midpoints=midpoints.psp(as.psp(network_lix))
  segment_lengths=lengths_psp(as.psp(network_lix))
  
  # Midpoints as a point pattern (on the original network)
  lpp_midpoints=lpp(midpoints,network)

  PAIs=matrix(NA,nrow=length(ks),ncol=length(ns))
  for (i in 1:length(ks)){
    # Compute distances between the middle points of those segments satisfying condition on k
    sig=which(rel_probs$probs>=mean(rel_probs$probs)+ks[i]*sd(rel_probs$probs))
    distances=crossdist.lpp(lpp_midpoints[sig,],X)
    for (j in 1:length(ns)){
      cat(paste0("k = ",ks[i],", n = ",ns[j]),"\n")
      hotspots=DRHotspots_k_n(X,rel_probs,ks[i],ns[j],event_distances=distances)
      if (!is.null(hotspots)){
        PAIs[i,j]=hotspots$PAI
      }
    }
  }
  rownames(PAIs)=paste0("k = ",ks)
  colnames(PAIs)=paste0("n = ",ns)
  return(PAIs)
}
