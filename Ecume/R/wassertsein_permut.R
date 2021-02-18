#' Permutation test based on Wasserstein distance
#'
#' @description Permutation test based on Wasserstein distance
#'
#' @param x Samples from the first distribution
#' @param y Samples from the second distribution. Only used if x is a vector.
#' @param iterations How many iterations to do to simulate the null distribution.
#' Default to 10^4.
#' @param fast If true, uses the \link[transport]{subwasserstein}
#' approximate function. Default to true if there are more than 1,000 samples
#' total.
#' @param S Number of samples to use in approximate mode. Must be set if `fast=TRUE`.
#' See \link[transport]{subwasserstein}.
#' @param ... Other parameters passed to \link[transport]{wasserstein} or
#' \link[transport]{wasserstein1d}
#' @examples
#'  x <- matrix(c(runif(100, 0, 1),
#'                runif(100, -1, 1)),
#'              ncol = 2)
#'  y <- matrix(c(runif(100, 0, 3),
#'                runif(100, -1, 1)),
#'              ncol = 2)
#'  # Set iterations to small number for runtime
#'  # Increase for more accurate results
#'  wasserstein_permut(x, y, iterations = 10^2)
#' @md
#' @return
#' A list containing the following components:
#' \itemize{
#'   \item *statistic* the Wasserstein distance between x and y.
#'   \item *p.value* the p-value of the permutation test.
#' }
#' @export
#' @importFrom dplyr if_else
#' @importFrom pbapply pblapply
#' @importFrom transport wasserstein wasserstein1d pp subwasserstein
#' @export
wasserstein_permut <- function(x, y, iterations = 10^4,
                               fast = nrow(x) + nrow(y) > 10^3,
                               S = NULL, ...) {
  if (ncol(x) != ncol(y)) stop("X and Y must have the same dimension")
  d <- ncol(x)
  nx <- nrow(x)
  ny <- nrow(y)
  X <- rbind(x, y)
  args <- list("a" = NULL, "b" = NULL, ...)
  if (d == 1) {
    dist_func <- transport::wasserstein1d
    trans <- identity
    fast <- FALSE
  } else {
    dist_func <- transport::wasserstein
    trans <- transport::pp
  }
  if (fast) {
    dist_func <- transport::subwasserstein
    args <- list("source" = NULL, "target" = NULL, ...)
    trans <- transport::pp
    args$S <- S
  }
  if (fast & !(is.numeric(S) && S > 0)) {
    stop("If fast mode is chosen, the S argument must be specified.")
  }
  args[[1]] <- trans(x); args[[2]] <- trans(y)
  og <- do.call(dist_func, args)
  null <- pbapply::pblapply(rep(0, iterations), function(rep) {
    permut <- sample(nx + ny)
    new_x <- X[permut[seq_len(nx)], ]
    new_y <- X[permut[seq(nx + 1, nx + ny)], ]
    args[[1]] <- trans(new_x); args[[2]] <- trans(new_y)
    return(do.call(dist_func, args))
  }) %>% unlist()
  pval <- max(1 / iterations, mean(og <= null))
  return(list("statistic" = og, "p.value" = pval))
}
