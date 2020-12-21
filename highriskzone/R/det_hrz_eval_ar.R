#' Determination of high-risk zone on smaller area of interest (evaluation area) than observation area.
#' 
#' \code{det_hrz_eval_ar} determines intensity based highriskzones if bomb crater observations are available
#' for a bigger area than the area of main interest (evaluation area).
#' All observations are used for intensity estimation, the highriskzone is however constructed only in the
#' evaluation area. Either based on specifying a failure probability alpha that indicates the probability of 
#' unobserved bombs outside the highriskzone but inside the evaluation area of interest (and not in the 
#' overall observation area) (criterion = "indirect"), or by specifying the threshold (maximum intensity of non-
#' exploded bombs outside the) highriskzone directly and intersecting the resulting hrz with the 
#' evaluation area (criterion = "direct").
#' @param ppdata  Observed spatial point process of class ppp in the observation area.
#' @param eval_ar area of interest specified via an object of class owin 
#' @param criterion  criterion to limit the high-risk zone, can be \code{"indirect"} (failure probability 
#' alpha) or \code{"direct"} (threshold, i.e. maximum intensity of unexploded bombs outside hrz) 
#' @param cutoff  Value of criterion (alpha or threshold)
#' @param intens  (optional) estimated intensity of the observed process (object of class "im") in (bigger) 
#' observation area, if not given, it will be estimated using \code{\link[spatstat.core]{density.ppp}}.
#' @param nxprob  Probability of having unobserved events.
#'                Default value is 0.1.
#' @param covmatrix  (optional) Covariance matrix of the kernel of a normal distribution, only needed for 
#'                    \code{type="intens"} if no intensity is given. If not given, it will be estimated
#'                    using \code{\link[ks]{Hscv}}. 
#' @return An object of class "\code{highriskzone}"
#' @export  
#' @examples
#' set.seed(12412)
#' spatstat.geom::spatstat.options(npixel=300)
#' data(craterB)
#' # reduce number of observations for faster computation
#' thin.craterB <- craterB[sample(1:craterB$n, 40)]
#' # define evaluation area of interest
#' eval.ar <- spatstat.geom::owin(xrange = c(0, 1900), yrange = c(0, 3400), 
#'                poly = matrix(c(250,250, 1200,1000,250,1000), byrow = TRUE, ncol = 2))
#'
#' hrzi1 <- det_hrz_eval_ar(thin.craterB, eval_ar = eval.ar, criterion = "direct",
#'                         cutoff = 3e-6, nxprob = .2)
#'
#' plot(hrzi1)
#' plot(thin.craterB, add = TRUE)
#' plot(eval.ar, add = TRUE)
#' plot(craterB$window, add = TRUE)

det_hrz_eval_ar<- function (ppdata, eval_ar, criterion = c("indirect", "direct"), cutoff,
                            intens = NULL, nxprob = 0.1, 
                            covmatrix = NULL) {
  win <- ppdata$window
  criterion <- match.arg(criterion, choices=c("indirect", "direct"))
  if (criterion=="indirect") {
    stopifnot(cutoff>0 & cutoff<1)
  } else {
    stopifnot(cutoff > 0)
  }
  stopifnot(nxprob>0 & nxprob<1)
  
  eval_ar <- as.owin(eval_ar)
  # Estimate Intensity
  if (is.null(intens)) {
    estim <- est_intens(ppdata, covmatrix = covmatrix)
    intens <- estim$intensest
    covmatrix <- estim$covmatrix
  }
  if (criterion == "indirect") {
  threshold <- det_threshold_eval_ar(intens = intens,
                                     eval_ar = eval_ar,
                                     alpha = cutoff, 
                                     nxprob = nxprob)
  } else {
    threshold <- (1-nxprob)/nxprob * cutoff
  }
  # Eval HRZ
  intens[setminus.owin(win, eval_ar)] <- 0
  HRZimage <- eval.im(intens > threshold)
  
  Rwindow <- owin(xrange = win$xrange, yrange = win$yrange, 
                  mask = as.matrix(HRZimage))
  threshold <- nxprob / (1-nxprob) * threshold
  result <- list(typehrz = "intens", criterion = paste(criterion, "_eval_ar", sep = ""), cutoff = cutoff, 
                 nxprob = nxprob, zone = Rwindow, threshold = threshold, calccutoff = NA, 
                 covmatrix = covmatrix)
  class(result) <- "highriskzone"
  return(result)
  
}

#' Determination of failure probability within evaluation area
#' @param intens estimated intensity
#' @param eval_ar evaluation area
#' @param threshold given threshold
#' @param nxprob constant probability of non-explosion   
det_alpha_eval_ar <- function(intens, eval_ar, threshold, nxprob = 0.1) {
  areaPixel <- intens$xstep * intens$ystep
  
  intens <- intens[as.owin(eval_ar), drop = FALSE]
  
  intensmatrix <- as.matrix(intens)
  intensmatrix[is.na(intensmatrix)] <- 0
  WwR <- (intensmatrix <= threshold)
  intens.intWwR <- sum(intensmatrix[WwR == 1] * areaPixel)
  alpha <- 1 - exp(-(nxprob/(1 - nxprob) * intens.intWwR))
  return(alpha)
}
#' Determination of necessary threshold to keep alpha in evaluation area
#' @param intens estimated intensity
#' @param eval_ar evaluation area
#' @param alpha desired failure probability in eval area
#' @param nxprob constant probability of non-explosion

det_threshold_eval_ar <- function(intens, eval_ar, alpha = 1e-05, nxprob = 0.1) {
  f <- function(logthreshold) {
    det_alpha_eval_ar(intens = intens, eval_ar = eval_ar, 
                      threshold = exp(logthreshold), nxprob = nxprob) - alpha
  }
  thres <- uniroot(f, lower = -100, upper = log(max(range(intens))))
  threshold <- exp(thres$root)
  return(threshold)
}

