#' Plots objects  produced by the \code{\link{logrr}} 
#' function.
#' 
#' Plots objects  of class \code{logrrenv} produced by the 
#' \code{\link{logrr}} function.
#' 
#' @param x An object of class \code{logrrenv}.
#' @param ... Additional graphical parameters passed to the 
#'   \code{\link[spatstat.geom]{image.im}} function.  See
#'   Details.
#' @param conlist Additional argument passed to the 
#'   \code{\link[spatstat.geom]{contour.im}} function.
#' @param main A main title for the plot.  Default is blank.
#' @details An important aspect of this plot is the
#'   color argument (\code{col}) used for displaying
#'   the regions outside the non-rejection envelopes.  If NULL
#'   (the implicit default), then the default color palette
#'   used by \code{\link[spatstat.geom]{image.im}} will be used. 
#'   Simpler schemes, e.g., c("blue", "white", "orange") can
#'   suffice. See the examples.
#' @method plot logrrenv
#' @seealso \code{\link[spatstat.geom]{plot.im}},
#'   \code{\link[spatstat.geom]{contour.im}}
#' @export
#' @examples
#' data(grave)
#' logrrsim = logrr(grave, nsim = 9)
#' plot(logrrsim)
#' # no border or ribben (legend).  Simple color scheme.
#' plot(logrrsim, col = c("blue", "white", "orange"), ribbon = FALSE, box = FALSE) 
#' # alternate color scheme
#' plot(logrrsim, col = topo.colors(12))
plot.logrrenv = function(x, ..., conlist = list(), main = "") {
  # if there were no simulations
  if (is.null(x$nrenv)) {
    spatstat.geom::image.im(x, ...)
  } else {
    # create temporary im object for plotting
    xtemp = spatstat.geom::im(mat = x$v, xcol = x$xcol, yrow = x$yrow)
    # determine which locations within non-rejection intervals or are NA
    which_na = rbind(which(x$nrenv$v == 0, arr.ind = TRUE),
                     which(is.na(x$v), arr.ind = TRUE))

    # arguments for contour.im function
    argc = list(x = xtemp, add = TRUE, conlist)
    
    # NA any locations outside window or inside non-rejection envelopes
    xtemp$v[which_na] = NA

    # arguments for image.im function
    argi = list(x = xtemp, main = main, ...)
       
    # plot colors for regions outside non-rejection envelopes
    do.call(spatstat.geom::image.im, argi)
    # add contour to plot
    do.call(spatstat.geom::contour.im, argc) 
  }
}
