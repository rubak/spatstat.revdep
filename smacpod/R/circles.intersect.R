#' Determine whether circles intersect
#' 
#' \code{circles.intersect} determines whether circles
#' intersect with each other.
#' 
#' The algorithm is based on the premise that two circles
#' intersect if, and only if, the distance between their
#' centroids is between the sum and the difference of their
#' radii.  I have squared the respective parts of the
#' inequality in the implemented algorithm.
#' 
#' @param coords A matrix of coordinates with the
#'   centroid of each circle.
#' @param r A vector containing the radii of the circles. 
#'   The length of \code{r} must equal the number of rows of 
#'   \code{coords}.
#' @return Returns a matrix of logical values indicating
#'   whether the circles intersect.
#' @author Joshua French
#' @export
#' @examples 
#' # first two circles intersect each other, 
#' # the next two circles intersect each other 
#' # (but not the previous ones)
#' # the last circles doesn't intersect any other circle
#' co = cbind(c(1, 2, 5, 6, 9), c(1, 2, 5, 6, 9))
#' r = c(1.25, 1.25, 1.25, 1.25, 1.25)
#' # draw circles
#' circles.plot(co, r)
#' # confirm intersections
#' circles.intersect(co, r)
#' 
#' # nested circles (don't intersect)
#' co = matrix(rep(0, 4), nrow = 2)
#' r = c(1, 1.5)
#' circles.plot(co, r)
#' circles.intersect(co, r)
circles.intersect <- function(coords, r) {
  # d = SpatialTools::dist1(coords)
  d = as.matrix(stats::dist(coords))
  if (length(r) != nrow(coords)) {
    stop("length(r) must be equal to nrow(coords)")
  }
  rmat1 = matrix(r, nrow = length(r), ncol = length(r))
  rmat2 = t(rmat1)
  lb = (rmat1 - rmat2)^2
  ub = (rmat1 + rmat2)^2
  return(lb <= d^2 & d^2 <= ub)
}

