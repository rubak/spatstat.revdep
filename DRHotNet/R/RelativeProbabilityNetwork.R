#' Computes the relative probability of observing a type of event along a linear network 
#' 
#' Given a marked point pattern lying on a linear network structure, this function uses kernel density estimation (KDE) to estimate a relative probability of occurrence for a type of event specified by the user through the marks of the pattern. The marks of a point pattern represent additional information of the events that are part of the pattern  
#' 
#' @param X - A \code{lpp} object representing a marked point pattern lying on a linear network (\code{linnet} object)
#' @param lixel_length - A numeric value representing a lixel length that will be used for creating a split version of the network contained in \code{X}. Then, the length of all the segments of the split network is below \code{lixel_length}
#' @param sigma - A numeric value representing the bandwidth parameter (in meters)
#' @param mark - Mark of \code{X} that is used to characterize the type of event. The algorithm searches microzones of the network where this mark is over- or underrepresented
#' @param category_mark - A numeric/character value from the set allowed in the chosen \code{mark} to compute the relative probability in relation to it
#' @return Returns a list that contains the relative probability values estimated along the network for the type of event specified by \code{mark} and \code{category_mark}
#' @examples 
#' library(DRHotNet)
#' library(spatstat.core)
#' library(spdep)
#' library(raster)
#' library(maptools)
#' \donttest{
#' rel_probs_rear_end <- RelativeProbabilityNetwork(X = SampleMarkedPattern, 
#' lixel_length = 50, sigma = 100, mark = "Collision", category_mark = "Rear-end")
#' }
#' @references Baddeley, A., Rubak, E., & Turner, R. (2015). Spatial point patterns: methodology and applications with R. Chapman and Hall/CRC.
#' @references Briz-Redon, A., Martinez-Ruiz, F., & Montes, F. (2019). Identification of differential risk hotspots for collision and vehicle type in a directed linear network. Accident Analysis & Prevention, 132, 105278.
#' @references Diggle, P. J. (2013). Statistical analysis of spatial and spatio-temporal point patterns. Chapman and Hall/CRC.
#' @references Kelsall, J. E., & Diggle, P. J. (1995). Kernel estimation of relative risk. Bernoulli, 1(1-2), 3-16.
#' @references McSwiggan, G., Baddeley, A., & Nair, G. (2017). Kernel density estimation on a linear network. Scandinavian Journal of Statistics, 44(2), 324-345.
#' @export
RelativeProbabilityNetwork <- function(X, lixel_length, sigma, mark, category_mark) {
  
  # Find position (column index) in marks(X) for the mark of interest
  
  marksX=data.frame(marks(X))
  if (!is.null(mark)){
    if (!is.null(names(marks(X)))){
      find_mark=which((names(marks(X))==mark)==T)
    } else{
      find_mark=1
    }
  }
  
  # Compute densities (equal-continuous, PDE method) as a function on a linear network
  
  # Considering all events
  density_function_all=as.linfun.linim(density.lpp(X, sigma = sigma))
  
  # Considering only the events whose mark is category_mark
  density_function_type=as.linfun.linim(density.lpp(X[marksX[,find_mark]==category_mark], sigma = sigma))
  
  # Extract network
  network=X$domain
  
  # Lixellize network
  if (lixel_length!=F){
    network_lix=lixellate(X$domain,eps=lixel_length)
    midpoints=midpoints.psp(as.psp(network_lix))
  } else{
    midpoints=midpoints.psp(as.psp(network))
  }
  
  # Midpoints as a point pattern (on the original network)
  lpp_midpoints=lpp(midpoints,network)
  
  # KDE computation
  density_values_all=density_function_all(lpp_midpoints$data$x,lpp_midpoints$data$y,
                     lpp_midpoints$data$seg,lpp_midpoints$data$tp)
  density_values_type=density_function_type(lpp_midpoints$data$x,lpp_midpoints$data$y,
                                  lpp_midpoints$data$seg,lpp_midpoints$data$tp)
  
  # output creation
  out=list()
  out$probs=density_values_type/density_values_all
  out$lixel_length=lixel_length
  out$sigma=sigma
  out$mark=mark
  out$category_mark=category_mark
  return(out)
}
