#' Plots an object obtained with \code{DiffHotspots_n_k}
#' 
#' This function plots a set of differential risk hotspots located along a linear network. An extension of the hotspots (including the kth order neighbourhs of the segments of the hotspots) is also plotted
#' 
#' @param X - A \code{lpp} object representing a marked point pattern lying on a road network (\code{linnet} object)
#' @param hotspots - A set of differential risk hotspots obtained with the function \code{DiffHotspots_n_k}
#' @param order_extension - A natural number indicating a neighbourhood order to be used for constructing an extension of the differential risk hotspots. The summary is also given for the segments forming this extension 
#' @param which.plot - A numeric vector indicating which differential risk hotspots to plot (according to the way they are ordered in \code{hotspots})
#' @param rotation_angle - A rotation angle (in degrees, from 0 to 180) to apply to the network (to improve visualization, if required). By default it is set to 0
#' @param eps_image - If set to \code{TRUE}, an .eps image is generated. By default it is set to \code{FALSE}
#' @examples 
#' library(DRHotNet)
#' library(spatstat)
#' library(spdep)
#' library(raster)
#' library(maptools)
#' \donttest{
#' rel_probs_rear_end <- RelativeProbabilityNetwork(X = SampleMarkedPattern, 
#' lixel_length = 50, sigma = 100, mark = "Collision", category_mark = "Rear-end")
#' hotspots_rear_end <- DRHotspots_k_n(X = SampleMarkedPattern, rel_probs = rel_probs_rear_end, 
#' k = 1, n = 30)
#' PlotHotspots(X = SampleMarkedPattern, hotspots = hotspots_rear_end)
#' }
#' @export
PlotHotspots <- function(X, hotspots, order_extension = NULL, which.plot = NULL, rotation_angle = 0, eps_image=F){
  
  network=X$domain
  lixel_length=hotspots$lixel_length
  sigma=hotspots$sigma
  mark=hotspots$mark
  category_mark=hotspots$category_mark
  k=hotspots$k
  n=hotspots$n
  
  if (is.null(order_extension)){
    order_extension=round(sigma/lixel_length)
  }
  
  if (hotspots$lixel_length!=F){
    network=lixellate(network,eps=hotspots$lixel_length)
    # project into the lixellized network
    X_aux=lpp(cbind(X$data$x,X$data$y),network)
    marks(X_aux)=marks(X)
    X=X_aux
  }
  network_lix=X$domain
  
  # Create sp object
  network_lix_sp=as.SpatialLines.psp(as.psp(X))
  
  # Extract hotspots segments
  
  # Filter if specified
  
  if (!is.null(which.plot)){
    segments_hotspots <- c()
    for (j in which.plot){
      segments_hotspots <- c(segments_hotspots, hotspots[[1]][[j]])
    }
  } else {
    segments_hotspots <- c()
    for (j in 1:length(hotspots[[1]])){
      segments_hotspots <- c(segments_hotspots, hotspots[[1]][[j]])
    }
  }
  
  # Neighbourhood matrix and hotspots extension
  
  W=NeighbourhoodMatrixNetwork(network_lix)
  segments_hotspots_extension=KthOrderNeighbours(segments_hotspots,W,order=order_extension)
  
  # Plot
  
  if (eps_image){
    setEPS()
    postscript(paste0("diff_risk_hotspots_k_",gsub("\\.","_",toString(k)),"_n_",n,"_lixel_",lixel_length,
                      "_sigma_",sigma,"_type_",mark,"_",category_mark,".eps"), family="Helvetica")
    par(mar=c(0,0,0,0))
    plot(elide(network_lix_sp, rotate=rotation_angle, center=apply(bbox(network_lix_sp), 1, mean)), col="black", lwd=1,
         main=paste0("Differential risk hotspots '",category_mark, "'",
                     " (", mark,"),", 
                     "\nlixel_length = ",lixel_length,", sigma = ",sigma, ",",
                     "\nk = ",k,", n = ",n), line=-4)
    plot(elide(network_lix_sp[segments_hotspots_extension,], rotate=rotation_angle, center=apply(bbox(network_lix_sp), 1, mean)),
         add=T,col="#fc9272",lwd=3)
    plot(elide(network_lix_sp[segments_hotspots,], rotate=rotation_angle, center=apply(bbox(network_lix_sp), 1, mean)),
         add=T,col="#de2d26",lwd=3)
    legend("bottom", legend=c("DRHotspot (center)",paste0("DRHotspot (extension), order = ", order_extension)),
           col=c("#de2d26","#fc9272"), lty=1, lwd=3, cex=0.85, title = NULL, horiz = T)
    dev.off()
  } else{
    par(xpd=TRUE)
    plot(elide(network_lix_sp, rotate=rotation_angle, center=apply(bbox(network_lix_sp), 1, mean)), col="black", lwd=1,
         main=paste0("Differential risk hotspots '",category_mark, "'",
                     " (", mark,"),", 
                     "\nlixel_length = ",lixel_length,", sigma = ",sigma, ",",
                     "\nk = ",k,", n = ",n))
    plot(elide(network_lix_sp[segments_hotspots_extension,], rotate=rotation_angle, center=apply(bbox(network_lix_sp), 1, mean)),
         add=T,col="#fc9272",lwd=3)
    plot(elide(network_lix_sp[segments_hotspots,], rotate=rotation_angle, center=apply(bbox(network_lix_sp), 1, mean)),
         add=T,col="#de2d26",lwd=3)
    legend("bottom", legend=c("DRHotspot (center)",paste0("DRHotspot (extension), order = ", order_extension)),
           col=c("#de2d26","#fc9272"), lty=1, lwd=3, cex=0.5, title = NULL, horiz = T)
  }
  

}