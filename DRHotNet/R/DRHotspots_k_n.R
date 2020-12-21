#' Identifies differential risk hotspots along a linear network given a vector of relative probabilities computed over the middle points of the segments of the network
#' 
#' Given a relative probability surface corresponding to the occurrence of a type of event along a linear network, this function filters and groups in hotspots those segments satisfying two conditions: 1) the relative probability in the segment exceeds the average relative probability per segment in \code{k} times the standard deviation of the complete set of probabilities estimated across all the segments of the network, and 2) there are \code{n} or more events at a distance below \code{sigma} from the middle point of the segment (\code{sigma} is obtained from the object \code{rel_probs} computed with the function \code{RelativeProbabilityNetwork}). In summary, \code{k} and \code{n} control the formation of differential risk hotspots along the network, given a set of relative probabilities covering the network. The choice of a higher value for \code{k} or \code{n} (or both) represents a more strict criterion and leads to a lower number of differential risk hotspots being detected. Users should test several values of \code{k} and \code{n} (sensitivity analysis on \code{k} and \code{n}) in order to reach reasonable choices for the research or practical purposes of their data analyses. This sensitivity analysis can be carried out with the \code{Sensitivity_k_n} function
#' 
#' @param X - A \code{lpp} object representing a marked point pattern lying on a road network (\code{linnet} object)
#' @param rel_probs - An object containing the relative probabilities of a specific type of event along the linear network contained in \code{X}, generated through the function \code{RelativeProbabilityNetwork}
#' @param k - A numeric value that controls the procedure of detecting differential risk hotspots (departure from average relative probability), as described above
#' @param n - A numeric value that controls the procedure of detecting differential risk hotspots (minimum size for the sample of events implicated in the computation of the relative probabilities), as described above
#' @param dist - A character indicating which distance to use. Two values are allowed: \code{path} (shortest-path distance) and \code{euclidean} (Euclidean distance). By default, the shortest-path distance is used. Change to \code{euclidean} to reduce the computation time or skip memory issues
#' @param event_distances - A matrix that contains the distances between the middle points of the segments satisfying the condition on parameter \code{k} and the events o \code{X}. By default it is set to \code{NULL}
#' @return Returns a list that contains the differential risk hotspots found for \code{X} and the type of event specified by \code{rel_probs}
#' @examples 
#' library(DRHotNet)
#' library(spatstat.core)
#' library(spdep)
#' library(raster)
#' library(maptools)
#' \donttest{
#' rel_probs_rear_end <- RelativeProbabilityNetwork(X = SampleMarkedPattern, 
#' lixel_length = 50, sigma = 100, mark = "Collision", category_mark = "Rear-end")
#' hotspots_rear_end <- DRHotspots_k_n(X = SampleMarkedPattern, rel_probs = rel_probs_rear_end, 
#' k = 1, n = 30)
#' }
#' @references Briz-Redon, A., Martinez-Ruiz, F., & Montes, F. (2019). Identification of differential risk hotspots for collision and vehicle type in a directed linear network. Accident Analysis & Prevention, 132, 105278.
#' @export
DRHotspots_k_n <- function(X,rel_probs,k,n,dist="path",event_distances=NULL){
  
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
  segment_lengths=lengths.psp(as.psp(network_lix))
  
  # Midpoints as a point pattern (on the original network)
  lpp_midpoints=lpp(midpoints,network)
  
  # sig initially contains the segments showing a high relative probability
  
  sig=which(rel_probs$probs>=mean(rel_probs$probs)+k*sd(rel_probs$probs))
  
  # find events close to the midpoints of sig, to determine n
  if (is.null(event_distances)){
    if (dist=="path"){
      distances=crossdist.lpp(lpp_midpoints[sig,],X)
    } else{
      X_planar=as.ppp(X)
      lpp_midpoints_planar=as.ppp(lpp_midpoints)
      distances=crossdist.ppp(lpp_midpoints_planar[sig,],X_planar)
    }
  } else{
    if (length(event_distances)==(length(sig)*length(X$data$x))){
      distances=event_distances
    } else{
      cat("The input event_distances has wrong dimensions\n")
      stop()
    }
  }
  closer_than_sigma=matrix(as.numeric(distances<rel_probs$sigma),nrow=nrow(distances))
  rownames(closer_than_sigma)=sig
  # print(sum(closer_than_sigma))
  
  # sig finally contains those segments that have n or more events at distance lower than sigma
  sig=sig[which(apply(closer_than_sigma,1,sum)>=n)]
  closer_than_sigma=closer_than_sigma[which(apply(closer_than_sigma,1,sum)>=n),]

  if (length(sig)>0){
    W=NeighbourhoodMatrixNetwork(network_lix)
    ordered_segments=c(1:length(sig))
    ordered_segments=c(1:length(sig))
    names(ordered_segments)=sig
    hotspot=rep(0,length(sig))
    count_zeros=length(which((hotspot==0)==T))
    current_hotspot=1
    while (count_zeros>0){
      ### new hotspot contains the first segment that has not been assigned to a hotspot yet
      start=which((hotspot==0)==T)[1]
      hotspot[start]=current_hotspot
      ### neighbouring segments are searched
      aux=as.numeric(names(ordered_segments)[start])
      aux_old=c()
      while (length(aux)>0){
        neigh=c()
        for (j in c(1:length(aux))){
          neigh=c(neigh,W$neighbours[[aux[j]]])
        }
        neigh=unique(neigh)
        if (length(aux_old>0)){
          find_old=c()
          for (j in c(1:length(aux_old))){
            find_old=c(find_old,which((neigh==aux_old[j])==T))
          }
          if (length(find_old)){
            neigh=neigh[-find_old]
          }
        }
        if (length(aux)>0){
          find=c()
          for (j in c(1:length(aux))){
            find=c(find,which((neigh==aux[j])==T))
          }
          if (length(find)>0){
            neigh=neigh[-find]
          }
        }
        ### update aux_old
        aux_old=aux
        ### keep neighbouring segments that should be part of the hotspot (those also in sig)
        aux=intersect(neigh,sig)
        ### assign hotspot number
        for (m in c(1:length(aux))){
          hotspot[ordered_segments[toString(aux[m])]]=as.numeric(current_hotspot)
        }
      }
      ### increase hotspot number
      current_hotspot=current_hotspot+1
      ### find segments not assigned yet
      count_zeros=length(which((hotspot==0)==T))
    }
    ### rearrange result
    result=list()
    segments_hotspots=c()
    for (j in c(1:length(unique(hotspot)))){
      result[[j]]=as.numeric(names(ordered_segments)[hotspot==j])
      segments_hotspots=c(segments_hotspots,as.numeric(names(ordered_segments)[hotspot==j]))
    }
    result_final=list()
    result_final[["DRHotspots"]]=result
    result_final$k=k
    result_final$n=n
    result_final$lixel_length=rel_probs$lixel_length
    result_final$sigma=rel_probs$sigma
    result_final$mark=rel_probs$mark
    result_final$category_mark=rel_probs$category_mark
    ### compute a global PAI
    marksX=as.data.frame(marks(X))
    if (!is.null(names(marks(X)))){
      index_mark=which(colnames(marksX)==result_final$mark)
    } else{
      index_mark=1
    }
    PAI=(length(which(marksX[X$data$seg%in%segments_hotspots,index_mark]==result_final$category_mark))/length(which(marksX[,index_mark]==result_final$category_mark)))/
      (sum(segment_lengths[segments_hotspots])/sum(segment_lengths))
    result_final$PAI_type=round(PAI,2)
  } else{
    result_final=NULL
    message("No differential risk hotspots found for the parameters provided. Tune k or n (or both) and try again")
  }
  return(result_final)
}
