#' Difference of estimated K functions
#' 
#' \code{kd} determines the difference in estimated K 
#' functions for a set of cases and controls.
#' 
#' This function relies internally on the 
#' \code{\link[spatstat.core]{Kest}} and 
#' \code{\link[spatstat.core]{eval.fv}} functions from the 
#' \code{spatstat} package.  The arguments are essentially
#' the same as the \code{\link[spatstat.core]{Kest}} function, 
#' and the user is referred there for more details about
#' the various arguments.
#' 
#' @param x A \code{ppp} object from the \code{spatstat} 
#'   package with marks for the case and control groups.
#' @param case The position of the name of the "case" group 
#'   in \code{levels(x$marks)}.  The default is 2.  
#'   \code{x$marks} is assumed to be a factor.  Automatic
#'   conversion is attempted if it is not.
#' @param domain Optional. Calculations will be restricted 
#'   to this subset of the window. See Details of 
#'   \code{\link[spatstat.core]{Kest}}.
#' @inheritParams spatstat.core::Kest
#'   
#' @return Returns an \code{fv} object.  See documentation 
#'   for \code{spatstat.core::Kest}.
#' @author Joshua French
#' @seealso \code{\link[spatstat.core]{Kest}}, 
#'   \code{\link[spatstat.core]{eval.fv}}
#' @references Waller, L.A. and Gotway, C.A. (2005). Applied
#'   Spatial Statistics for Public Health Data. Hoboken, NJ:
#'   Wiley.
#' @export
#' @examples 
#' data(grave)
#' kd = kd(grave)
#' plot(kd)
kd = function(x, case = 2, r = NULL, rmax = NULL, 
              breaks = NULL, 
              correction = c("border", "isotropic", "Ripley", "translate"), 
              nlarge = 3000, domain = NULL, 
              var.approx = FALSE, ratio = FALSE) {
  cases = which(x$marks == levels(x$marks)[case])
  K_case = spatstat.core::Kest(x[cases, ], r = r, rmax = rmax, breaks = breaks, correction = correction, nlarge = nlarge,
                Domain = domain, var.approx = var.approx, ratio = ratio)
  K_control = spatstat.core::Kest(x[-cases, ], r = r, rmax = rmax, breaks = breaks, correction = correction, nlarge = nlarge,
                   domain = domain, var.approx = var.approx, ratio = ratio)
  spatstat.core::eval.fv(K_case - K_control)
}
