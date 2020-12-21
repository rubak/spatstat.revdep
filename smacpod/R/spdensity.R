#' Kernel smoothed spatial density of point pattern
#' 
#' \code{spdensity} computes a kernel smoothed spatial
#' density function from a point pattern.  This function is
#' basically a wrapper for \code{\link[spatstat.core]{density.ppp}}.
#' The \code{\link[spatstat.core]{density.ppp}} function computes
#' the spatial intensity of a point pattern; the \code{spdensity}
#' function scales the intensity to produce a true spatial density. 
#' @inheritParams spatstat::density.ppp
#' 
#' @return This function produces the spatial density of \code{x}
#' as an object of class \code{im} from the \code{spatstat} package.  
#' @author Joshua French
#' @export
#' @seealso \code{\link[spatstat.core]{density.ppp}}
#' @references Waller, L.A. and Gotway, C.A. (2005).  Applied Spatial Statistics for Public Health Data.  Hoboken, NJ: Wiley.
#' @examples 
#' data(grave)
#' contour(spdensity(grave))

spdensity <- function(x, sigma = NULL, ..., weights=NULL, edge = TRUE, 
                     varcov = NULL, at = "pixels", 
                     leaveoneout = TRUE, adjust = 1, diggle = FALSE,
                    kernel = "gaussian", scalekernel = is.character(kernel),
                     positive = FALSE, verbose = TRUE) {
  d <- spatstat.core::density.ppp(x = x, sigma = sigma, ..., weights = weights,
              edge = edge, varcov = varcov, at = at, leaveoneout = leaveoneout,
              adjust = adjust, diggle = diggle,
              se = FALSE, kernel = kernel, scalekernel = scalekernel,
              positive = positive, verbose = verbose)
  d$const <- spatstat.geom::integral.im(d)
  d$v <- d$v/d$const
  class(d) <- c(class(d), "spdensity")
  return(d)
}
