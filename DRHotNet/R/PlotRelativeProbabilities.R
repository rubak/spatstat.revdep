#' Plots an object obtained with \code{RelativeProbabilityNetwork}
#' 
#' This function plots the relative probability of occurrence of a type of event along a linear network
#' 
#' @param X - A \code{lpp} object representing a marked point pattern lying on a road network (\code{linnet} object)
#' @param rel_probs - An object containing the relative probabilities of a specific type of event along the linear network contained in \code{X}, generated through the function \code{RelativeProbabilityNetwork}
#' @param rotation_angle - A rotation angle (in degrees, from 0 to 180) to apply to the network (to improve visualization, if required). By default it is set to 0
#' @param eps_image - If set to \code{TRUE}, an .eps image is generated. By default it is set to \code{FALSE}
#' @examples 
#' library(DRHotNet)
#' library(spatstat)
#' library(spdep)
#' library(maptools)
#' \donttest{
#' rel_probs_rear_end <- RelativeProbabilityNetwork(X = SampleMarkedPattern, 
#' lixel_length = 50, sigma = 100, mark = "Collision", category_mark = "Rear-end")
#' PlotRelativeProbabilities(X = SampleMarkedPattern, rel_probs = rel_probs_rear_end)
#' }
#' @export
PlotRelativeProbabilities <- function(X, rel_probs, rotation_angle = 0, eps_image=F){
  
  floor_dec=function(x, level=1) round(x - 5*10^(-level-1), level)
  ceiling_dec=function(x, level=1) round(x + 5*10^(-level-1), level)
  
  network=X$domain
  lixel_length=rel_probs$lixel_length
  sigma=rel_probs$sigma
  mark=rel_probs$mark
  category_mark=rel_probs$category_mark
  
  if (rel_probs$lixel_length!=F){
    network=lixellate(network,eps=rel_probs$lixel_length)
    # project into the lixellized network
    X_aux=lpp(cbind(X$data$x,X$data$y),network)
    marks(X_aux)=marks(X)
    X=X_aux
  }
  network_lix=X$domain
  
  # Create sp object
  network_lix_sp=as.SpatialLines.psp(as.psp(X))
  
  # Plot
  
  if (eps_image){
    setEPS()
    postscript(paste0("rel_probs_lixel_",lixel_length,"_sigma_",sigma,"_type_",mark,"_",category_mark,".eps"), family="Helvetica")
    par(mar=c(0,0,0,0))
    breaks_colors=c(-1,quantile(rel_probs$probs,0.2),quantile(rel_probs$probs,0.4),
                    quantile(rel_probs$probs,0.6),quantile(rel_probs$probs,0.8),1)
    color=cut(rel_probs$probs,breaks = breaks_colors)
    palette_edit=c("#2C7BB6","#ABD9E9","#bababa","#FDAE61","#D7191C")
    palette(palette_edit)
    plot(elide(network_lix_sp, rotate=rotation_angle, center=apply(bbox(network_lix_sp), 1, mean)),col=color,lwd=2,
         main=paste0("Relative probabilities '",category_mark, "'",
                     " (", mark,"),", 
                     "\nlixel_length = ",lixel_length,", sigma = ",sigma), line=-2) 
    
    legend("bottom", legend=c(paste0("[",floor_dec(min(rel_probs$probs),2),", ",round(breaks_colors[2],2),"]"),
                              paste0("]",round(breaks_colors[2],2),", ",round(breaks_colors[3],2),"]"),
                              paste0("]",round(breaks_colors[3],2),", ",round(breaks_colors[4],2),"]"), 
                              paste0("]",round(breaks_colors[4],2),", ",round(breaks_colors[5],2),"]"),
                              paste0("]",round(breaks_colors[5],2),", ",ceiling_dec(max(rel_probs$probs),2),"]")),
           col=palette_edit, lty=1, lwd=3, cex=0.85, title = NULL, horiz = T)
    dev.off()
  } else{
    par(xpd=TRUE)
    breaks_colors=c(-1,quantile(rel_probs$probs,0.2),quantile(rel_probs$probs,0.4),
                    quantile(rel_probs$probs,0.6),quantile(rel_probs$probs,0.8),1)
    color=cut(rel_probs$probs,breaks = breaks_colors)
    palette_edit=c("#2C7BB6","#ABD9E9","#bababa","#FDAE61","#D7191C")
    palette(palette_edit)
    plot(elide(network_lix_sp, rotate=rotation_angle, center=apply(bbox(network_lix_sp), 1, mean)),col=color,lwd=2,
         main=paste0("Relative probabilities '",category_mark, "'",
                     " (", mark,"),", 
                     "\nlixel_length = ",lixel_length,", sigma = ",sigma)) 
    
    legend("bottom", legend=c(paste0("[",floor_dec(min(rel_probs$probs),2),", ",round(breaks_colors[2],2),"]"),
                              paste0("]",round(breaks_colors[2],2),", ",round(breaks_colors[3],2),"]"),
                              paste0("]",round(breaks_colors[3],2),", ",round(breaks_colors[4],2),"]"), 
                              paste0("]",round(breaks_colors[4],2),", ",round(breaks_colors[5],2),"]"),
                              paste0("]",round(breaks_colors[5],2),", ",ceiling_dec(max(rel_probs$probs),2),"]")),
           col=palette_edit, lty=1, lwd=3, cex=0.5, title = NULL, horiz = T)
  }
  
  
}