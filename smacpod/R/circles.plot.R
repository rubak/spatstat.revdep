#' Plot circles
#' 
#' \code{plot.circles} creates a plot with one or more 
#' circles (or adds them to an existing plot).
#' 
#' @param coords A matrix of coordinates with the centroid
#'   of each circle.
#' @param r A vector containing the radii of the circles. 
#'   The length of \code{r} must equal the number of rows of
#'   \code{coords}.
#' @param add A logical value indicating whether the circles
#'   should be added to an existing plot.  Default is
#'   \code{FALSE}.
#' @param ... Additional arguments passed to the 
#'   \code{\link[graphics]{plot}} function.
#' @inheritParams plotrix::draw.circle
#' @param border A vector with the desired border of each
#'   circle.  The length should either be 1 (in which case
#'   the border is repeated for all circles) or should match
#'   the number of rows of \code{coords}.
#' @param ccol A vector with the desired color of each
#'   circle.  The length should either be 1 (in which case
#'   the color is repeated for all circles) or should match
#'   the number of rows of \code{coords}.
#' @param clty A vector with the desired line type of each
#'   circle.  The length should either be 1 (in which case
#'   the line type is repeated for all circles) or should match
#'   the number of rows of \code{coords}.
#' @param density A vector with the density for a patterned fill.
#'   The length should either be 1 (in which case
#'   the density is repeated for all circles) or should match
#'   the number of rows of \code{coords}. See \code{\link[graphics]{polygon}}
#' @param angle A vector with the angle of a patterned fill.
#'   The length should either be 1 (in which case
#'   the angle is repeated for all circles) or should match
#'   the number of rows of \code{coords}. See \code{\link[graphics]{polygon}}
#' @param clwd A vector with the desired line width of each
#'   circle.  The length should either be 1 (in which case
#'   the line width is repeated for all circles) or should match
#'   the number of rows of \code{coords}.
#' @return NULL
#' @author Joshua French
#' @importFrom graphics plot
#' @importFrom plotrix draw.circle
#' @seealso \code{\link[plotrix]{draw.circle}}, \code{\link[graphics]{polygon}}
#' @export
#' @examples 
#' co = cbind(c(1, 2, 5, 6, 9), c(1, 2, 5, 6, 9))
#' r = c(1.25, 1.25, 1.25, 1.25, 1.25)
#' # draw circles
#' circles.plot(co, r)
#' circles.plot(co, r, 
#'    ccol = c("blue", "blue", "orange", "orange", "brown"),
#'    density = c(10, 20, 30, 40, 50),
#'    angle = c(45, 135, 45, 136, 90))
circles.plot <- function(coords, r, add = FALSE, ..., 
                         nv = 100, border = NULL,
                         ccol = NA,
                         clty = 1,
                         density = NULL,
                         angle = 45,
                         clwd = 1) {
  dots <- list(...)
  if (is.null(dim(coords))) {
    stop("coords must be a matrix-like object")
  }
  if (length(r) != nrow(coords)) {
    stop("length(r) must be equal to nrow(coords)")
  }
  if (length(add) != 1) {
    stop("add should be a single logical value")
  }
  if (!is.logical(add)) {
    stop("add should be a single logical value")
  }
  if (length(nv) == 1) {
    nv <- rep(nv, length(r))
  } else if (length(nv) != length(r)) {
    stop("nv must have length 1 or length(r)")
  }
  if (is.null(border)) {
    border <- vector("list", length(r))
  } else if (length(border) == 1) {
    border <- rep(border, length(r))
  } else if (length(border) != length(r)) {
    stop("border must have length 1 or length(r)")
  } else if (length(border) == length(r)) {
    border <- as.list(border)
  }
  if (length(ccol) == 1) {
    ccol <- rep(ccol, length(r))
  } else if (length(ccol) != length(r)) {
    stop("ccol must have length 1 or length(r)")
  }
  if (length(clty) == 1) {
    clty <- rep(clty, length(r))
  } else if (length(clty) != length(r)) {
    stop("clty must have length 1 or length(r)")
  }
  if (is.null(density)) {
    density <- vector("list", length(r))
  } else if (length(density) == 1) {
    density <- as.list(rep(density, length(r)))
  } else if (length(density) != length(r)) {
    stop("density must have length 1 or length(r)")
  } else if (length(density) == length(r)) {
    density <- as.list(density)
  }
  if (length(angle) == 1) {
    angle <- rep(angle, length(r))
  } else if (length(angle) != length(r)) {
    stop("angle must have length 1 or length(r)")
  }
  if (length(clwd) == 1) {
    clwd <- rep(clwd, length(r))
  } else if (length(clwd) != length(r)) {
    stop("clwd must have length 1 or length(r)")
  }
  
  if (!add) {
    # determine bounds of plot
    lowx <- min(coords[,1] - r)
    highx <- max(coords[,1] + r)
    lowy <- min(coords[,2] - r)
    highy <- max(coords[,2] + r)
    dots$x <- c(lowx, highx)
    dots$y <- c(lowy, highy)
    dots$type <- "n"
    if (is.null(dots$xlim)) {
      dots$xlim <- c(lowx, highx)
    }
    if (is.null(dots$ylim)) {
      dots$ylim <- c(lowy, highy)
    }
    if (is.null(dots$asp)) {
      dots$asp <- 1
    }
    if (is.null(dots$xlab)) {
      dots$xlab <- "x"
    }
    if (is.null(dots$ylab)) {
      dots$ylab <- "y"
    }
    do.call(plot, dots)
  }
  for (i in seq_along(r)) {
    plotrix::draw.circle(coords[i, 1], coords[i, 2], r[i],
                         nv = nv[i], border = border[[i]], col = ccol[i],
                         lty = clty[i], density = density[[i]],
                         angle = angle[i], lwd = clwd[i])
  }
}

