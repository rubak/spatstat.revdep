#' Performs a summary of a set of differential risk hotspots located along a linear network
#' 
#' This function provides a basic summary of each differential risk hotspot provided in the object \code{hotspots} passed to the function. This includes the proportion of the type of event in each hotspot, the total length of the hotspot, a type-specific prediction accuracy index (\code{PAI_type}). Furthermore, this summary is also provided for an extension of each of the hotspots
#' 
#' @param X - A \code{lpp} object representing a marked point pattern lying on a linear network (\code{linnet} object)
#' @param rel_probs - An object containing the relative probabilities of a specific type of event along the linear network contained in \code{X}, generated through the function \code{RelativeProbabilityNetwork}
#' @param hotspots - A set of differential risk hotspots obtained with the function \code{DiffHotspots_n_k}
#' @param order_extension - A natural number indicating a neighbourhood order to be used for constructing an extension of the differential risk hotspots. The summary is also given for the segments forming this extension 
#' @param compute_p_value - A Boolean value allowing the user to compute a p-value representing the statistical significance of each differential risk hotspot 
#' @param n_it - Number of simulations performed for the estimation of the p-value (if \code{compute_p_value} = \code{T})
#' @return Returns a \code{data.frame} providing a summary of a set of differential risk hotspots. Each row of the output corresponds to one hotspot
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
#' hotspots_summary <- SummaryDRHotspots(X = SampleMarkedPattern, rel_probs = rel_probs_rear_end, 
#' hotspots = hotspots_rear_end)
#' }
#' @export
SummaryDRHotspots <- function(X,rel_probs,hotspots,order_extension=NULL,compute_p_value=F,n_it=40){
  
  network=X$domain
  sigma=rel_probs$sigma
  lixel_length=rel_probs$lixel_length
  mark=rel_probs$mark
  category_mark=rel_probs$category_mark
  
  if (is.null(order_extension)){
    order_extension=round(sigma/lixel_length)
  }
  
  if (rel_probs$lixel_length!=F){
    network=lixellate(network,eps=rel_probs$lixel_length)
    # project into the lixellized network
    X_aux=lpp(cbind(X$data$x,X$data$y),network)
    marks(X_aux)=marks(X)
    X=X_aux
  }
  network_lix=X$domain
  middle_points=midpoints.psp(as.psp(network_lix))
  segment_lengths=lengths_psp(as.psp(network_lix))
  
  # Neighbourhood matrix
  
  W=NeighbourhoodMatrixNetwork(network_lix)

  result=c()
  rownames_aux=rownames(hotspots[[2]])
  for (i in c(1:length(hotspots[[1]]))){
    segments_hotspots=hotspots[[1]][[i]]
    segments_hotspots_extension=KthOrderNeighbours(segments_hotspots,W,order=order_extension)
    if (compute_p_value==T){
      total=n_it
      pb <- winProgressBar(min = 0, max = total, label=sprintf(paste0("Differential risk hotspot: ",i)))
      rel_probs_sim=c()
      for (n in c(1:n_it)){
        setWinProgressBar(pb, n, label=sprintf(paste0("Differential risk hotspot: ",i)))
        ### Mark permutation
        X_perm=MarkPermutation(X)
        aux=RelativeProbabilityNetwork(X_perm,lixel_length=rel_probs$lixel_length,
                                       sigma=rel_probs$sigma,mark=mark,category_mark=category_mark)
        ### Compute average intensity in the hotspot 
        weighted_average_sim=sum(aux$probs[segments_hotspots]*(segment_lengths[segments_hotspots]/sum(segment_lengths[segments_hotspots])))
        rel_probs_sim=c(rel_probs_sim,weighted_average_sim)
      }
      rel_probs_sim=sort(rel_probs_sim)
      weighted_average_true=sum(rel_probs$probs[segments_hotspots]*(segment_lengths[segments_hotspots]/sum(segment_lengths[segments_hotspots])))
      p_value_est=1-length(which(rel_probs_sim<=weighted_average_true))/length(rel_probs_sim)
    }
    marksX=as.data.frame(marks(X))
    if (!is.null(names(marks(X)))){
      index_mark=which(colnames(marksX)==mark)
    } else{
      index_mark=1
    }
    PAI=(length(which(marksX[X$data$seg%in%segments_hotspots,index_mark]==category_mark))/length(which(marksX[,index_mark]==category_mark)))/
      (sum(segment_lengths[segments_hotspots])/sum(segment_lengths))
    PAI_extension=(length(which(marksX[X$data$seg%in%segments_hotspots_extension,index_mark]==category_mark))/
                   length(which(marksX[,index_mark]==category_mark)))/
                   (sum(segment_lengths[segments_hotspots_extension])/sum(segment_lengths))
    ### Save results
    row_result=c()
    row_result[1]=i
    row_result[2]=round(length(which(marksX[X$data$seg%in%segments_hotspots,index_mark]==category_mark)),2)
    row_result[3]=round(length(which(X$data$seg%in%segments_hotspots)),2)
    row_result[4]=round(row_result[2]/row_result[3],2)
    row_result[5]=round(sum(segment_lengths[segments_hotspots]),2)
    row_result[6]=round(PAI,2)
    row_result[7]=round(length(which(marksX[X$data$seg%in%segments_hotspots_extension,index_mark]==category_mark)),2)
    row_result[8]=round(length(which(X$data$seg%in%segments_hotspots_extension)),2)
    row_result[9]=round(row_result[7]/row_result[8],2)
    row_result[10]=round(sum(segment_lengths[segments_hotspots_extension]),2)
    row_result[11]=round(PAI_extension,2)
    if (compute_p_value==T){
      row_result[12]=round(p_value_est,6)
    }
    if (compute_p_value==T){
      names(row_result)=c("DRHotspot",
                          "Events type (ctr)","All events (ctr)",
                          "Prop. (ctr)","Length in m (ctr)","PAI_type (ctr)",
                          "Events type (ext)","All events (ext)",
                          "Prop. (ext)","Length in m (ext)","PAI_type (ext)",
                          "p-value")
    } else{
      names(row_result)=c("DRHotspot",
                          "Events type (ctr)","All events (ctr)",
                          "Prop. (ctr)","Length in m (ctr)","PAI_type (ctr)",
                          "Events type (ext)","All events (ext)",
                          "Prop. (ext)","Length in m (ext)","PAI_type (ext)")
    }
    result=rbind(result,row_result)
  }
  result=as.data.frame(result)
  rownames(result)=NULL
  return(result)
}
