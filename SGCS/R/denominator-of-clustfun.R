#' Estimate the denominator of clustering function
#' 
#'
#' @param x Point pattern
#' @param r Vector of distances to estimate the function
#' @param correction Border correction. Either "border", "none" or "best"(="border).
#' @param ... Ignored.
#' 
#' @details
#' Support in 3D therefore only for cuboidal windows. (inc. rotated)
#' 
#' @return
#' \code{\link{fv}}-object.
#' 
#' @useDynLib SGCS
#' @import spatstat
#' @export

clustfun_denominator <- function(x, r, correction="best", ...) {
  x <- internalise_pp(x)
  ### range
  r <- default_r(x, r)
  
  ### Distances for speed
  x$pairwise_distances <- pairwise_distances(x)
  
  ### Border distances for correction
  correction_i <- correction %in% c("border","best")
  x$edgeDistances <- if(correction_i) edge_distance(x) else rep(max(r), x$n)
  
  
  ### Compute:
  res <- .External("SGCS_clustfun_denominator_c",
                   x,
                   r,
                   PACKAGE="SGCS"
  )
  
  A <- if(x$dim==2) pi*r^2 else pi*r^3 * 4/3
  # 
  lam <- x$n/x$area
  lpr <- lam * A
  #p  <- (1 - exp(-lpr) * (lpr + 1)  )
  theo <- 0.5 * lam^2 * A^2
  
  # make fv suitable
  c.final<-fv( data.frame(r=r, theo=theo, cd=res),
               argu = "r",
               alim = range(r),
               ylab = substitute(c_denom(r), NULL),
               desc = c("distance argument r", "Theoretical values unknown", "Denom. of clustering function"),
               valu = "cd",
               fmla = ".~r",
               fname="cd"  
  )
  
  c.final
}
