#' Determine nearest neighbors
#' 
#' \code{nn} determines the nearest neighbors for a set of 
#' observations based on a distance matrix.
#' 
#' This function determine nearest neighbors in two ways: 1.
#' number of neighbors or 2. distance.
#' 
#' If \code{method = "c"}, then \code{k} specifies the total
#' number of neighbors to return for each observation.  
#' 
#' If \code{method = "d"}, then \code{k} specifies the maximum
#' distance for which an observation is considered a 
#' neighbor.  
#' 
#' The function returns the neighbors for each 
#' observation.
#' 
#' @param d A square distance matrix for the 
#'   coordinates of interest.
#' @param k The number of neighbors to return (if 
#'   \code{method = "c"}) or the distance for which 
#'   observations are considered neighbors (if \code{method 
#'   = "d"}).
#' @param method The method of determining the neighbors. 
#'   The default is \code{"c"}, specifying that the \code{k}
#'   nearest neighbors (the number of neighbors) for each 
#'   observation should be returned.  The alternative is 
#'   \code{"d"}, meaning that neighbors are determined by 
#'   their distance from an observation.  In that case, two 
#'   observations are neighbors if their separation distance
#'   is less or equal to \code{k}.
#' @param self A logical indicating whether an observation 
#'   is a neighbor with itself.  The default is 
#'   \code{FALSE}.
#'   
#' @return Returns a list with the nearest neighbors of each
#'   observation. For each element of the list, the indices
#'   order neighbors from nearest to farthest.
#' @author Joshua French
#' @export
#' @examples 
#' data(grave)
#' # make distance matrix
#' d = as.matrix(dist(cbind(grave$x, grave$y)))
#' # 3 nearest neighbors
#' nnc = nn(d, k = 3, method = "c")
#' # nearest neighbors within k units of each observation
#' nnd = nn(d, k = 200, method = "d")
nn <- function(d, k, method = "c", self = FALSE) {
  if (!is.matrix(d)) {
    stop("d must be a matrix of coordinates")
  }
  if (nrow(d) != ncol(d)) {
    stop("d should be a square matrix")
  }
  if (!is.numeric(k) || k <= 0 || length(k) > 1) {
    stop("k must be a single posive numeric value")
  }
  if (!is.element(method, c("c", "d"))) {
    stop("valid options for method are 'c' or 'd'")
  }
  if (!is.logical(self)) {
    stop("self must be a logical value")
  }
  nws = (self + 1) %% 2 # switch TRUE to 0, FALSE to 1
  if ((k + nws >= ncol(d)) & method == "c") {
    stop("k is larger than available number of neighbors")
  }
  # split d into list of columns of d
  d <- split(d, rep(1:ncol(d), each = nrow(d)))
  
  if (method == "c") {
    # implement slightly different behavior based on whether the observation 
    # is a neighbor to itself
    # nws was converted to numeric to shift index appropriately, depending
    # on what is wanted
    # ordered distances (column represent results for each observation)
    # select only needed number of rows
    # not that we increment the returned indexes returned by 1
    # if each observation is not a neighbor with itself
    lapply(d, function(x) order(x)[(1 + nws):(k + nws)])
  }else {
    # for each row in d: order the distances, only return the
    # ones with distances less than k.  Since distances are ordered, 
    # this will simply be the first q ordered elements,
    # where q is the number of observations in that row 
    # with distances <= k.  Adjustment made based on whether
    # whether we should include the observation itself (self = TRUE)
    lapply(d, function(x) order(x)[(1 + nws):(sum(x <= k))])
  }
}
