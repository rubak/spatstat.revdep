#' @title Calibration Simplex
#' @aliases CalSim
#'
#' @description Generates an object of class \code{calibration_simplex} which can be used to assess the calibration
#' of ternary probability forecasts. The Calibration Simplex can be seen as generalization of the reliability diagram
#' for binary probability forecasts. For details on the interpretation of the calibration simplex, see Wilks (2013). Be
#' aware that some minor changes have been made compared to the calibration simplex as suggested by Wilks (2013) (see note below).
#' 
#' As a somewhat experimental feature, multinomial p-values can be used for uncertainty quantification, that is, as a tool
#' to judge whether the observed discrepancies may be merely coincidental or whether the predictions may in fact be miscalibrated, see Resin (2020, Section 4.2).
#'
#' @param n A natural number.
#' @param p1 A vector containing the forecasted probabilities for the first (1) category, e.g. below-normal.
#' @param p2 A vector containing the forecasted probabilities for the second (2) category, e.g. near-normal.
#' @param p3 A vector containing the forecasted probabilities for the third (3) category, e.g. above-normal.
#' @param obs A vector containing the observed outcomes (Categories are encoded as 1 (e.g. below-normal), 2 (e.g. near-normal) and 3 (e.g. above-normal)).
#' @param test_stat A string indicating which test statistic is to be used for the multinomial test in each bin. 
#' Options are "LLR" (log-likelihood ratio; default), "Chisq" (Pearson's chi-square) and "Prob" (probability mass statistic). See details
#' @param percentagewise Logical, specifying whether probabilities are percentagewise (summing to 100) or not (summing to 1).
#' 
#' @return A list with class "calibration_simplex" containing
#'   \item{\code{n}}{As input by user or default.}
#'   \item{\code{n_bins}}{Computed from \code{n}. Number of hexagons.}
#'   \item{\code{n_obs}}{Total number of observations.}
#'   \item{\code{freq}}{Vector of length \code{n_bins} containing the number of observations within each bin.}
#'   \item{\code{cond_rel_freq}}{Matrix containing the observed outcome frequencies within each bin.}
#'   \item{\code{cond_ave_prob}}{Matrix containing the average forecast probabilities within each bin.}
#'   \item{\code{pvals}}{Exact multinomial p-values within each bin. See details.}
#'
#' @rdname calibration_simplex
#' @export
#'
#' @details Only two of the three forecast probability vectors (\code{p1}, \code{p2} and \code{p3}) need to be specified.
#' 
#' The p-values are based on multinomial tests comparing the observed frequencies within a bin 
#' with the average forecast probabilities within the bin as outlined in Resin (2020, Section 4.2).
#' The p-values are exact and do not rely on asymptotics, however, it is assumed that the true 
#' distribution (under the hypothesis of forecast calibration) within each bin 
#' is approximated well by the multinomial distribution. If \code{n} is small the 
#' approximation may be poor, resulting in unreliable p-values. p-Values less than 0.0001 are not
#' exact but merely indicate a value less than 0.0001.
#'
#' @examples
#' attach(ternary_forecast_example)   #see also documentation of sample data
#' #?ternary_forecast_example
#'
#' # Calibrated forecast sample
#' calsim0 = calibration_simplex(p1 = p1, p3 = p3, obs = obs0)
#' plot(calsim0,use_pvals = TRUE) # with multinomial p-values
#'
#' # Overconfident forecast sample
#' calsim1 = calibration_simplex(p1 = p1, p3 = p3, obs = obs1)
#' plot(calsim1)
#'
#' # Underconfident forecast sample
#' calsim2 = calibration_simplex(p1 = p1, p3 = p3, obs = obs2)
#' plot(calsim2,use_pvals = TRUE) # with multinomial p-values
#'
#' # Unconditionally biased forecast sample
#' calsim3 = calibration_simplex(p1 = p1, p3 = p3, obs = obs3)
#' plot(calsim3)
#'
#' # Using a different number of bins
#' calsim = calibration_simplex(n=4, p1 = p1, p3 = p3, obs = obs3)
#' plot(calsim)
#'
#' calsim = calibration_simplex(n=13, p1 = p1, p3 = p3, obs = obs3)
#' plot(calsim,               # using some additional plotting parameters:
#'      error_scale = 0.5,    # errors are less pronounced (smaller shifts)
#'      min_bin_freq = 100,   # dots are plotted only for bins,
#'                            # which contain at least 100 forecast-outcome pairs
#'      category_labels = c("below-normal","near-normal","above-normal"),
#'      main = "Sample calibration simplex")
#'
#' detach(ternary_forecast_example)


calibration_simplex = function(n,p1,p2,p3,obs,test_stat,percentagewise){
  UseMethod("calibration_simplex")
}
