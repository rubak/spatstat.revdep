#' Log ratio of spatial densities
#' 
#' \code{logrr} computes the log ratio of spatial density 
#' functions for cases and controls.  The numerator in this 
#' ratio is related to the "cases" and the denominator to 
#' the "controls".  If \code{nsim > 0}, then pointwise 
#' non-rejection envelopes are used to assess potential 
#' clustering of cases and controls relative to each other.
#' 
#' The \code{plot} function makes it easy to visualize the
#' log ratio of spatial densities (if \code{nsim = 0}) or
#' the regions where the log ratio deviates farther from
#' than what is expected under the random labeling 
#' hypothesis (i.e., the locations of potential clustering).
#' The shaded regions indicate the locations of potential
#' clustering.
#' 
#' @inheritParams spatstat.core::density.ppp
#' @param sigma Standard deviation of isotropic smoothing
#'   kernel for cases. Either a numerical value, or a function that
#'   computes an appropriate value of \code{sigma}.
#' @param sigmacon Standard deviation of isotropic smoothing
#'   kernel for controls.  Default is the same as 
#'   \code{sigma}.
#' @param case The position of the name of the "case" group 
#'   in \code{levels(x$marks)}.  The default is 2.  
#'   \code{x$marks} is assumed to be a factor.  Automatic
#'   conversion is attempted if it is not.
#' @param nsim The number of simulated data sets from which 
#'   to construct the non-rejection envelopes under the random
#'   labeling hypothesis.  Default is 0 (i.e., no
#'   envelopes).
#' @param level The level used for the pointwise
#'   non-rejection envelopes.
#' @param alternative The direction of the significance test
#'   to identify potential clusters using a Monte Carlo test
#'   based on the pointwise non-rejection envelopes.  Default is
#'   \code{"two.sided"} (logrr != 0).  The values \code{"less"}
#'   (logrr < 0) and \code{"greater"} (logrr > 0) are also valid.
#' @param bwargs A list of arguments for the bandwidth 
#'   function supplied to \code{sigma} and \code{sigmacon},
#'   if applicable.
#'   
#' @return The function produces an object of type 
#'   \code{logrrenv}.  Its components are similar to those 
#'   returned by the \code{density.ppp} function from the 
#'   \code{spatstat} package, with the intensity values 
#'   replaced by the log ratio of spatial densities of f and
#'   g.  Includes an array \code{simr} of dimension c(nx, 
#'   ny, nsim + 1), where nx and ny are the number of x and 
#'   y grid points used to estimate the spatial density. 
#'   \code{simr[,,1]} is the log ratio of spatial densities 
#'   for the observed data, and the remaining \code{nsim} 
#'   elements in the third dimension of the array are the 
#'   log ratios of spatial densities from a new ppp 
#'   simulated under the random labeling hypothesis.
#' @details The \code{two.sided} alternative test assesses 
#'   whether the observed ratio of log densities deviates 
#'   more than what is expected under the random labeling 
#'   hypothesis.  When the test is significant, this 
#'   suggests that the cases and controls are clustered, 
#'   relative to the other.  The \code{greater} alternative 
#'   assesses whehter the cases are more clustered than
#'   the controls.  The \code{less} alternative
#'   assesses whether the controls are more clustered than 
#'   the cases.  If the estimated density of the case or 
#'   control group becomes too small, this function may 
#'   produce warnings due to numerical underflow. Increasing
#'   the bandwidth (sigma) may help.
#' @author Joshua French (and a small chunk by the authors 
#'   of the \code{\link[spatstat.core]{density.ppp}}) function 
#'   for consistency with the default behavior of that 
#'   function)
#' @export
#' @references Waller, L.A. and Gotway, C.A. (2005). Applied
#'   Spatial Statistics for Public Health Data. Hoboken, NJ:
#'   Wiley.  
#'   
#'   Kelsall, Julia E., and Peter J. Diggle. "Kernel
#'   estimation of relative risk." Bernoulli (1995): 3-16. 
#'   
#'   Kelsall, Julia E., and Peter J. Diggle. "Non-parametric
#'   estimation of spatial variation in relative risk." 
#'   Statistics in Medicine 14.21-22 (1995): 2335-2342.
#' @examples 
#' data(grave)
#' r = logrr(grave)
#' plot(r)
#' r2 = logrr(grave, sigma = spatstat.core::bw.scott)
#' plot(r2)
#' rsim = logrr(grave, nsim = 9)
#' plot(rsim)
logrr = function(x, sigma = NULL, sigmacon = NULL, case = 2, 
                 nsim = 0, level = 0.90, alternative = "two.sided", ..., 
                 bwargs = list(), weights = NULL, edge = TRUE, 
                 varcov = NULL, at = "pixels", leaveoneout = TRUE, 
                 adjust = 1, diggle = FALSE, 
                 kernel = "gaussian",
                 scalekernel = is.character(kernel),
                 positive = FALSE, verbose = TRUE) {
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
  if (level <= 0 | level >= 1) stop("level must be between 0 and 1.")
  if (!is.element(alternative, c("two.sided", "greater", "lower"))) stop("alternative is not valid.")
  alpha = 1 - level

  # determine bandwidth if necessary
  if (is.function(sigma)) { # use user-supplied function, if given
    which_bwargs <- which(names(bwargs) %in% names(formals(sigma))[-1])
    if (length(which_bwargs) > 0 ) {
      sigma = do.call(sigma, c(list(X = x, bwargs[which_bwargs])))
    } else {
      sigma = do.call(sigma, list(X = x))
    }
  }
  if (is.null(sigma)) { # use spatstat.core::bw.relrisk if nothing given
    which_bwargs <- names(bwargs) %in% names(formals(spatstat.core::bw.relrisk))[-1]
    if (length(which_bwargs) > 0 ) {
      sigma = do.call(spatstat.core::bw.relrisk, c(list(X = x, bwargs[which_bwargs])))
    } else {
      sigma = do.call(spatstat.core::bw.relrisk, list(X = x))
    }
  }
  if (is.null(sigmacon)) sigmacon = sigma # determine sigmacon, if NULL
  
  cases = which(x$marks == levels(x$marks)[case])
  N1 = length(cases)
  r = spdensity(x = x[cases,], sigma = sigma, ..., 
                weights = weights[cases],
                edge = edge, varcov = varcov, at = at, 
                leaveoneout = leaveoneout,
                adjust = adjust, diggle = diggle,
                kernel = kernel, 
                scalekernel = scalekernel,
                positive = positive, verbose = verbose)
  
  g = spdensity(x = x[-cases,], sigma = sigmacon, ..., 
                weights = weights[-cases],
                edge = edge, varcov = varcov, at = at, 
                leaveoneout = leaveoneout,
                adjust = adjust, diggle = diggle,
                kernel = kernel, 
                scalekernel = scalekernel,
                positive = positive, verbose = verbose)
  r$v = log(r$v) - log(g$v)
  r$nrenv = NULL
  
  if (nsim > 0) {
    simr2 <- pbapply::pblapply(seq_len(nsim), function(i) {
      cases = sample(x$n, N1)
      fsim = spdensity(x = x[cases,], sigma = sigma, ..., 
                       weights = weights[cases],
                       edge = edge, varcov = varcov, at = at, 
                       leaveoneout = leaveoneout,
                       adjust = adjust, diggle = diggle,
                       kernel = kernel, 
                       scalekernel = scalekernel,
                       positive = positive, verbose = verbose)
      
      gsim = spdensity(x = x[-cases,], sigma = sigmacon, ..., 
                       weights = weights[-cases],
                       edge = edge, varcov = varcov, at = at, 
                       leaveoneout = leaveoneout,
                       adjust = adjust, diggle = diggle,
                       kernel = kernel, 
                       scalekernel = scalekernel,
                       positive = positive, verbose = verbose)
      
      log(fsim$v) - log(gsim$v)
    })
    simr2[[nsim + 1]] <- simr2[[1]]
    simr2[[1]] <- r$v
    simr2 <- abind::abind(simr2, along = 3)
    
    r$simr = simr2
    r$nrenv = nrenv(r, level = level, alternative = alternative)
    class(r) = c("logrrenv", class(r))
  }
  
  r$window = x$window
  return(r)
}


