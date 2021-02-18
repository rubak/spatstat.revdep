.full_kernel <- function(x, y, norms, kernel, ...) {
  z <- rbind(x, y)
  K <- kernlab::kernelFast(kernel, z, z, norms)
  return(K)
}

MMD2u <- function(K, m, n) {
  # The MMD^2_u unbiased statistic.
  Kx2 <- K[seq_len(m), seq_len(m)]
  diag(Kx2) <- NA
  Ky2 <- K[(m + 1):(m + n), (m + 1):(m + n)]
  diag(Ky2) <- NA
  Kxy2 <- K[seq_len(m), (m + 1):(m + n)]
  term1 <- (1.0 / (m * (m - 1.0))) * sum(Kx2, na.rm = TRUE)
  term2 <- (1.0 / (n * (n - 1.0))) * sum(Ky2, na.rm = TRUE)
  term3 <- (2.0 / (m * n)) * sum(Kxy2)
  return(term1 + term2 - term3)
}

compute_null_distribution_u <- function(K, m, n, iterations = 10000) {
  # Compute the bootstrap null-distribution of MMD2u.
  mmd2u_null <- pbapply::pblapply(seq_len(iterations), function(i) {
    idx <- sample(m + n)
    K_i <- K[idx, idx]
    return(MMD2u(K_i, m, n))
  })
  return(mmd2u_null)
}

.ind_kernels <- function(x, y, m, n, kernel, frac = 1){
  # Sample for X
  l <- min(round(sqrt(frac * m)), round(sqrt(frac * n)))
  idx <- sample(m, 2 * l)
  x1 <- x[idx[seq_len(l)], ]
  x1 <- as.matrix(x1)
  norms_x1 <- rowSums(x1^2)
  x2 <- x[idx[seq(l + 1, 2 * l)], ]
  x2 <- as.matrix(x2)
  norms_x2 <- rowSums(x2^2)
  # Sample for Y
  l <- round(sqrt(frac * n))
  idy <- sample(n, 2 * l)
  y1 <- y[idy[seq_len(l)], ]
  y1 <- as.matrix(y1)
  norms_y1 <- rowSums(y1^2)
  y2 <- y[idx[seq(l + 1, 2 * l)], ]
  y2 <- as.matrix(y2)
  norms_y2 <- rowSums(y2^2)
  # Kernels
  Kx <- kernlab::kernelFast(kernel, x1, x2, norms_x1) %>% as.vector()
  Ky <- kernlab::kernelFast(kernel, y1, y2, norms_y1) %>% as.vector()
  Kxy <- kernlab::kernelFast(kernel, x1, y1, norms_x1) %>% as.vector()
  Kyx <- kernlab::kernelFast(kernel, y2, x2, norms_y2) %>% as.vector()
  Kx_ <- c(Kx, Kxy)
  Ky_ <- c(Ky, Kyx)
  return(list("Kx_" = Kx_, "Ky_" = Ky_, "l" = l))
}

MMDl <- function(Kx_, Ky_,  l) {
  # The MMD^2_u linear statistic.
  term1 <- (1 / l) * sum(Kx_[seq_len(l)])
  term2 <- (1 / l) * sum(Ky_[seq_len(l)])
  term3 <- (1 / l) * sum(Kx_[seq(l + 1, 2* l)])
  term4 <- (1 / l) * sum(Ky_[seq(l + 1, 2* l)])
  return(term1 + term2 - term3 - term4)
}

compute_null_distribution_l <- function(sample_Ks, iterations = 10000) {
  # Compute the bootstrap null-distribution of MMDl.
  l <- sample_Ks$l
  Kx_ <- sample_Ks$Kx_
  Ky_ <- sample_Ks$Ky_
  mmdl_null <- pbapply::pblapply(seq_len(iterations), function(i){
    ids <- sample(2 * l)
    Kx_i <- Kx_[ids]
    Ky_i <- Ky_[ids]
    return(MMDl(Kx_i, Ky_i, l))
  })
  return(unlist(mmdl_null))
}

#' Perform the Maximum Mean Discrepancy unbiased bootstrap test
#'
#' @description Maximum Mean Discrepancy Unbiased Test
#'
#' @param x d-dimensional samples from the first distribution
#' @param y d-dimensional samples from the first distribution
#' @param kernel A character that must match a known kernel. See details.
#' @param type Which statistic to use. One of 'unbiased' or 'linear'. See
#' Gretton et al for details. Default to 'unbiased' if the two vectors are of
#' length less than \code{1000} and to 'linear' otherwise.
#' @param null How to asses the null distribution. This can only be set to exact
#' if the `type` is 'unbiased' and the `kernel` is 'rbf'.
#' @param iterations How many iterations to do to simulate the null distribution.
#' Default to 10^4. Only used if `null` is 'permutations'
#' @param frac For the linear statistic, how many points to sample. See details.
#' @param ... Further arguments passed to kernel functions
#' @examples
#' x <- matrix(rnorm(1000, 0, 1), ncol = 10)
#' y <- matrix(rnorm(1000, 0, 2), ncol = 10)
#' mmd_test(x, y)
#' mmd_test(x, y, type = "linear")
#' x <- matrix(rnorm(1000, 0, 1), ncol = 10)
#' y <- matrix(rnorm(1000, 0, 1), ncol = 10)
#'  # Set iterations to small number for runtime
#'  # Increase for more accurate results
#' mmd_test(x, y, iterations = 10^2)
#' @details
#' This computes the MMD^2u unbiased statistic or the MMDl linear statistic
#' from Gretton et al. The code relies on the pairwise_kernel function from the
#' python module sklearn. To list the available kernels, see the examples.
#' @md
#' @references
#' Gretton, A., Borgwardt, K., Rasch, M. J., Schölkopf, B., & Smola, A. (2012).
#' *A Kernel Two-Sample Test* Journal of Machine Learning Research (2012)
#' @return
#' A list containing the following components:
#' \itemize{
#'   \item *statistic* the value of the test statistic.
#'   \item *p.value* the p-value of the test.
#' }
#' @importFrom pbapply pbapply
#' @import kernlab
#' @export
mmd_test <- function(x, y, kernel = 'rbfdot',
                     type = ifelse(min(nrow(x), nrow(y)) < 1000,
                                   "unbiased", "linear"),
                     null = c("permutation", "exact"),
                     iterations = 10^3,
                     frac = 1,
                     ...) {
  null <- match.arg(null)
  if (null == "exact" && (type == "linear" | kernel != 'rbfdot')) {
    stop("The exact mode only works with the unbiased statistic and the rbf kernel")
  }
  if (is.character(kernel)) {
    args <- list(...)
    kernel <- do.call(kernel, args)
  }
  # Compute MMD^2_u or MMDl, its null distribution and the p-value of the
  # kernel two-sample test.
  m <- nrow(x)
  n <- nrow(y)
  if (type == "unbiased") {
    norms <- c(rowSums(x^2), rowSums(y^2))
    K <- .full_kernel(x, y, norms, kernel, ...)
    statistic <- MMD2u(K, m, n)
    if (null == "permutation") {
      nulls <- compute_null_distribution_u(K, m, n, iterations = iterations)
      p.value <- max(1 / iterations, mean(nulls > statistic))
    }
  }
  if (type == "linear") {
    sample_Ks <- .ind_kernels(x, y, m , n, kernel, frac = frac, ...)
    statistic <- MMDl(sample_Ks$Kx_, sample_Ks$Ky_, sample_Ks$l)
    nulls <- compute_null_distribution_l(sample_Ks, iterations = iterations)
    p.value <- max(1 / iterations, mean(nulls > statistic))
  }
  return(list("statistic" = statistic, "p.value" = p.value))
}
