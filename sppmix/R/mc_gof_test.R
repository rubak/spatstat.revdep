#' Monte Carlo goodness of fit test
#'
#' @description
#' Performs a Monte Carlo test of goodness-of-fit for a given point pattern.
#' The entertained model is a Poisson with mixture of normals intensity
#' surface.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #mc_gof}
#'
#' @param pp Point pattern object of class \code{ppp}.
#' @param intsurf Object of class \code{intensity_surface}.
#' @param alpha Significance level for the goodness-of-fit test.
#' @param truncate Requests to truncate the components
#' of the mixture intensity to have all their mass
#' within the window of the intensity object intsurf. Default is FALSE.
#' @param L Number of iterations requested; default is 20000.
#' @param burnin Number of initial realizations to discard. By default, it is 1/10 of the total number of iterations.
#' @details The test statistic is the average of the average distances between the
#' points assigned to the jth mixture component from the mean of the component.
#' The Monte Carlo test utilizes realizations from the posterior
#' predictive distribution to obtain the critical point,
#' i.e., the \code{a}th percentile of the distribution of the test statistic.
#' Make sure that L is large in order to get accurate results.
#' @author Jiaxun Chen, Sakis Micheas
#' @seealso \code{\link{normmix}},
#' \code{\link{rsppmix}}
#' @examples
#' \donttest{
#' # Create the intensity surface
#' intsurf1 <- normmix(ps = c(.3, .7), mus = list(c(0.2, 0.2), c(.8, .8)), sigmas =
#'  list(.01*diag(2), .01*diag(2)), lambda = 100, win = spatstat.geom::square(1))
#' # Generate a point pattern
#' pp1 <- rsppmix(intsurf1)
#' # Assess goodness-of-fit. Since this is the right model, we should get gof. Make
#' # sure L is large for more accurate results
#' mc_gof(pp1, intsurf1, 0.05)
#' # Create another intensity surface
#' intsurf2 <- normmix(ps = c(.5, .5), mus = list(c(0.2, 0.8), c(.8, .2)), sigmas =
#'  list(.01*diag(2), .01*diag(2)), lambda = 100, win = spatstat.geom::square(1))
#' # Assess goodness-of-fit against this Poisson. Since this is the wrong model,
#' # we should NOT get gof
#' mc_gof(pp1, intsurf2, 0.05)}
#'
#' @export
mc_gof <- function(pp, intsurf, alpha=0.5, L = 20000, burnin = floor(.1*L),
                   truncate = FALSE) {
  # check first
  if (L < 1000) {
    stop("L needs to be larger than 1000.")
  }
  n <- pp$n
  if (n == 0) {
    stop("No points in the point pattern.")
  }
  m <- intsurf$m
  if (m == 0) {
    stop("No components in the mixture.")
  }
  zeronj <- 0
  win <- intsurf$window
  # get posterior realizations
  post <- est_mix_damcmc(pp = pp, m = m, truncate = truncate, L = L)
  T_mean <- rep(0, (L - burnin))
  # get posterior predictive sample
  for (i in 1:(L-burnin)) {
    lambda <- post$genlamdas[i+burnin]
    ps <- post$genps[(i+burnin), ]/sum(post$genps[(i+burnin), ])
    mus <- split(post$genmus[, , (i+burnin)], 1:m)
    sigma <- post$gensigmas[(i+burnin), ]
    mix_real <- normmix(ps, mus, sigma, lambda = lambda, win = win)
    ow <- options("warn")
    options(warn = -1)
    pp_pred <- rsppmix(mix_real, truncate = truncate)
    options(ow)
    # compute the test statistics for the predicted pattern
    summean <- 0
    ind <- unique(pp_pred$comp)
    if (length(ind) != m) {
      zeronj <- zeronj + 1
    }
      for (k in ind) {
        ppk <- pp_pred[pp_pred$comp == k]
        distk <- sqrt(rowSums((sweep(cbind(ppk$x, ppk$y), 2, mus[[k]]))^2))
        summean <- summean + mean(distk)
      } # end loop of k
    T_mean[i] <- summean/m
  } # end look of i
  if (zeronj != 0) {
    message(paste("There are", zeronj, "cases with components have 0 points in them.\nMC test could be reporting wrong.") )
  }
  sortT_mean <- sort(T_mean)
  # compute test statistcs for the given pattern
  summean <- 0
  probz <- GetAvgLabelsDiscrete2Multinomial_sppmix(
    post$genzs[(burnin + 1):L, ], m)
    for (j in 1:m) {
    distj <- probz[, j]*sqrt(rowSums((sweep(cbind(pp$x, pp$y), 2, intsurf$mus[[j]]))^2))
    summean <- summean + mean(distj)
  } # end loop of j
  T_meanTS <- summean/m
  # test for T_mean
  T_meanalpha <- sortT_mean[floor((1 - alpha)*(L - burnin))]
  # Ts_mean <- T_meanTS > T_meanalpha
  # if (Ts_mean == TRUE) {
  #   result_mean <- "One-sided GOF using T_meanTS, Reject Null.\\n
  #   Mixture Model DOES NOT FIT WELL"
  # } else {
  #   result_mean <- "One-sided GOF using T_meanTS, Cannot Reject Null.\\n
  #   Mixture Model FITS WELL"
  # }
  p_mean <- mean(T_meanTS < T_mean)
  if (p_mean < alpha) {
    result_mean <- "\nSmall p-value, Mixture Model does not fit well.\n "
  } else {
    result_mean <- "\nLarge p-value, Mixture Model fits the data well. \n"
  }
  result <- paste("\n         Monte Carlo Goodness of Fit Test \nTest Statistic:",
                  round(T_meanTS, 4), " Critical Value:", round(T_meanalpha, 4),
                  "\nNull Hypothesis: Mixture Model fits the pattern well.",
                  "\nAlternative: Mixture Model is not appropriate. \np-value:",
                  round(p_mean, 4), result_mean)
  cat(result)
}
