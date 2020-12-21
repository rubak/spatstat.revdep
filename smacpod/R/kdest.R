#' Difference of estimated K functions
#' 
#' \code{kdest} determines the difference in estimated K
#' functions for a set of cases and controls.  Non-rejection
#' envelopes can also be produced.
#' 
#' This function relies internally on the 
#' \code{\link[spatstat]{Kest}} and 
#' \code{\link[spatstat]{eval.fv}} functions from the 
#' \code{spatstat} package.  The arguments are essentially
#' the same as the \code{\link[spatstat]{Kest}} function, 
#' and the user is referred there for more details about
#' the various arguments.
#' 
#' @param x A \code{\link[spatstat]{ppp}} object 
#'   package with marks for the case and control groups.
#' @param case The position of the name of the "case" group 
#'   in \code{levels(x$marks)}.  The default is 2. 
#'   \code{x$marks} is assumed to be a factor.  Automatic 
#'   conversion is attempted if it is not.
#' @param nsim An non-negative integer.  Default is 0.  The
#'   difference in estimated K functions will be calculated
#'   for \code{nsim} data sets generated under the random
#'   labeling hypothesis.  These will be used to construct
#'   the non-rejection envelopes.
#' @param level The level used for the non-rejection envelopes. 
#'   Ignored if \code{nsim} is 0.
#' @param domain Optional. Calculations will be restricted
#'   to this subset of the window. See Details of
#'   \code{\link[spatstat]{Kest}}.
#' @inheritParams spatstat::Kest
#'   
#' @return Returns a \code{kdenv} object.  See documentation
#'   for \code{spatstat::Kest}.
#' @author Joshua French
#' @export
#' @seealso \code{\link[spatstat]{Kest}},
#'   \code{\link[spatstat]{eval.fv}}
#' @references Waller, L.A. and Gotway, C.A. (2005). 
#'   Applied Spatial Statistics for Public Health Data. 
#'   Hoboken, NJ: Wiley.
#' @examples 
#' data(grave)
#' kd1 = kdest(grave)
#' plot(kd1, iso ~ r, ylab = "difference", legend = FALSE, main = "")
#' kd2 = kdest(grave, nsim = 9, level = 0.8)
#' plot(kd2)
kdest = function(x, case = 2, nsim = 0, level = 0.95, r = NULL, 
                 rmax = NULL, breaks = NULL, 
                 correction = c("border", "isotropic", "Ripley", "translate"), 
                 nlarge = 3000, domain = NULL, 
                 var.approx = FALSE, ratio = FALSE) {
  if (!is.element("ppp", class(x))) stop("x must be a ppp object")
  if (is.null(x$marks)) stop("x must be marked as cases or controls")
  if (!is.factor(x$marks)) {
    message("converting marks(x) to a factor")
    x$marks <- factor(x$marks)
  }
  if (!is.factor(x$marks)) stop("The marks(x) must be a factor")
  nlev = length(levels(x$marks))
  if (case < 1 || case > nlev) stop("case must be an integer between 1 and length(levels(x$marks))")
  if (nsim < 0 | !is.finite(nsim)) stop("nsim must be a non-negative integer")
  if (length(level) != 1) stop("level must have length 1")
  if (level <= 0 | level >= 1) stop("level should be between 0 and 1")
  
  if (nsim == 0) {
    out = kd(x, case = case, 
             r = r, rmax = rmax, breaks = breaks, 
             correction = correction, 
             nlarge = nlarge, domain = domain, 
             var.approx = var.approx, ratio = ratio)
    out = list(out = out)
  } else {
    #min/max envelope
    out = spatstat::envelope(x, kd, case = case, nsim = nsim, 
                             savefuns = TRUE, 
                             simulate = expression(spatstat::rlabel(x, permute = TRUE)), 
                             r = r, rmax = rmax, 
                             breaks = breaks, 
                             correction = correction, 
                             nlarge = nlarge, 
                             domain = domain, 
                             var.approx = var.approx, 
                             ratio = ratio)
    simfuns <- as.data.frame(attr(out, "simfuns"))
    simfuns[,1] <- out$obs
    l = apply(simfuns, 1, stats::quantile, prob  = (1 - level)/2)
    u = apply(simfuns, 1, stats::quantile, prob = 1 - (1 - level)/2)
    out = list(out = out, qlo = l, qhi = u)
  }
  out$r = out$out$r
  out$case_label = levels(x$marks)[case]
  out$control_label = levels(x$marks)[-case]
  out$nsim = nsim
  out$level = level
  class(out) = "kdenv"
  return(out)
}