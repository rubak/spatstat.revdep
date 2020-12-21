#' Determine non-overlapping clusters
#' 
#' Determine non-overlapping clusters from a list of
#' potential clusters.
#' 
#' The function takes a list of potential clusters. Each
#' element of the list contains a potential cluster. The
#' potential clusters are defined by the location indices of
#' the regions comprising the clusters.  Starting with the
#' first potential cluster, the function excludes every
#' potential cluster that intersects the first (any
#' potential cluster that shares indices).  Moving onto the
#' next non-overlapping cluster, the process is repeated. 
#' The function returns the indices (in the list of
#' clusters) of the clusters that do not overlap.
#' 
#' @param x A list containing the potential clusters.
#' @return A vector with the list indices of the 
#'   non-overlapping clusters.
#' @author Joshua French
#' @export
#' @examples 
#' x = list(1:2, 1:3, 4:5, 4:6, 7:8)
#' noc(x)
noc = function(x) {
  if (!is.list(x)) {
    stop("x should be a list with location ids for the regions in the clusters")
  }
  remain_idx = seq_along(x)
  i = 1
  u = 1
  while (i < length(x)) {
    # do clusters intersect?
    inter = unlist(lapply(x[remain_idx], function(q) length(intersect(q, x[[i]])) == 0), use.names = FALSE)
    
    # non-intersecting clusters
    remain_idx = remain_idx[which(inter)]
    if (length(remain_idx) > 0) {
      i = min(remain_idx)
      u = c(u, i)
    } else {
      i = length(x) + 1
    }
  }
  return(u)
}
