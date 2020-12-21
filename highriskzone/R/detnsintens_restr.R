#' Determination of the intensity for the Neyman-Scott simulation.
#'
#' Used in function bootcor_restr.
#'
#' @param ppdata  observed point pattern whose estimated intensity (adjusted for
#'          thinning and divided by "clustering") is used for simulating the
#'          parent process
#' @param radius  radius of the circles around the parent points in which the cluster
#'          points are located
#' @param weights Vector of observation probabilities associated with the observations contained in \code{ppdata}.
#' @return A pixel image (object of class "im"). See \code{\link[spatstat.core]{density.ppp}}.
#' @seealso \code{\link[spatstat.core]{density.ppp}}, \code{\link[spatstat.geom]{boundingbox}},
#'          \code{\link[spatstat.geom]{owin}}, \code{\link[ks]{Hscv}}




det_nsintens_restr <- function(ppdata, radius, weights){

  pppsim <- ppdata
  frameWin <- boundingbox(pppsim$window)
  dilatedWin <- owin(frameWin$xrange + 1.2*c(-radius, radius), frameWin$yrange + 1.2*c(-radius, radius))
  pppsim$window <- dilatedWin

  covmatrixsim <- Hscv(cbind(pppsim$x, pppsim$y))
  intensest <- density(pppsim, varcov=covmatrixsim, weights=weights, positive=TRUE)

  return(intensest)
}
