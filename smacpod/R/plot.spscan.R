#' Plots object from \code{\link{spscan.test}}.
#' 
#' Plots object of class \code{scan} from
#' \code{\link{spscan.test}}.
#' 
#' If \code{border}, \code{ccol}, \code{clty}, or 
#' \code{clwd} are specified, then the length of these
#' vectors must match \code{nrow(x$coords)}.
#' 
#' @param x An object of class \code{spscan}.
#' @param ... Additional graphical parameters passed to the
#'   \code{\link[spatstat.geom]{plot.ppp}} function.
#' @param nv The number of verticies to draw the cluster 
#' circles. Default is 100.
#' @param border The border color of the circle.  Default is
#'   NULL, meaning black.
#' @param ccol Fill color of the circles.  Default is NULL,
#'   indicating empty.
#' @param clty Line type of circles.  Default is NULL,
#'   indicting \code{lty = 1}.
#' @param clwd Line width of circles.  Default is NULL,
#'   indicating \code{lwd = 2} for the most likely cluster and 
#'   \code{lwd = 1} for the rest.
#' @method plot spscan
#' @seealso \code{\link[spatstat.geom]{plot.ppp}}, \code{\link[plotrix]{draw.circle}}
#' @export
#' @examples
#' data(grave)
#' out = spscan.test(grave, case = 2, alpha = 0.1)
#' plot(out, chars = c(1, 20), main = "most likely cluster")
plot.spscan = function(x, ..., nv = 100, border = NULL, 
                     ccol = NULL, clty = NULL, clwd = NULL) {
  if (class(x) != "spscan") {
    stop("x should be a spscan object from spscan.test function")
  }
  spatstat.geom::plot.ppp(x$ppp, ...)
  
  # number of centroids
  nc = length(x$clusters)
  
  # set default values
  if (is.null(border)) border = rep(1, nc)
  if (is.null(clty)) clty = rep(1, nc)
  if (is.null(ccol)) ccol = rep(NA, nc)
  if (is.null(clwd)) clwd = c(2, rep(1, nc - 1))
  
  # more sanity checking
  if (length(border) != nc) stop("if specified, border must have length equal to nrow(x$coords)")
  if (length(clty) != nc) stop("if specified, border must have length equal to nrow(x$coords)")
  if (length(ccol) != nc) stop("if specified, border must have length equal to nrow(x$coords)")
  if (length(clwd) != nc) stop("if specified, border must have length equal to nrow(x$coords)")
  
  # plot clusters
  for (i in seq_len(nc)) {
    plotrix::draw.circle(x$clusters[[i]]$coords[1, 1], x$clusters[[i]]$coords[1, 2], 
                         x$clusters[[i]]$r, nv = nv, border = border, col = ccol[i], 
                         lty = clty[i], lwd = clwd[i])
  }
}
