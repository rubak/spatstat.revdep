.check_d <- function(x, y = NULL) {
  if (is.null(y)) {
    check <- dplyr::n_distinct(lapply(x, ncol) %>% unlist)
    check <- check == 1
  } else {
    check <- ncol(x) == ncol(y)
  }
  if (!check) {
    stop("All distributions must have the same dimensions")
  }
}

#' Classifier k-sample test
#'
#' @description Classifier k-sample test
#'
#' @param x Samples from the first distribution or a list of samples
#' from k distribution
#' @param y Samples from the second distribution. Only used if x is a vector.
#' @param split How to split the data between training and test. Default to .7
#' @param thresh Value to add to the null hypothesis. See details.
#' @param method Which model(s) to use during training. Default to knn.
#' @param control Control parameters when fitting the methods.
#' See \link[caret]{trainControl}
#' @param ... Other parameters passed to \link[caret]{train}
#' @examples
#'  x <- matrix(c(runif(100, 0, 1),
#'                runif(100, -1, 1)),
#'              ncol = 2)
#'  y <- matrix(c(runif(100, 0, 3),
#'                runif(100, -1, 1)),
#'              ncol = 2)
#'  classifier_test(x, y)
#' @details
#' See Lopez-Paz et .al for more background on those tests.
#' @md
#' @references
#' Lopez-Paz, D., & Oquab, M. (2016). Revisiting Classifier Two-Sample Tests, 1–15. Retrieved from http://arxiv.org/abs/1610.06545
#' @return
#' A list containing the following components:
#' \itemize{
#'   \item *statistic* the value of the test statistic.
#'   \item *p.value* the p-value of the test.
#' }
#' @export
#' @importFrom dplyr bind_rows n_distinct group_by sample_n ungroup
#' @importFrom methods is
#' @import caret e1071
#' @importFrom stats pbinom
classifier_test <- function(x, y, split = .7, thresh = 0,
                            method = "knn",
                            control = caret::trainControl(method = "cv"),
                            ...) {
  if ("list" %in% methods::is(x)) {
    .check_d(x)
    names(x) <- paste0("C", seq_along(x))
    x <- lapply(x, as.data.frame)
    X <- dplyr::bind_rows(x, .id = "type")
    X$type <- as.factor(X$type)
  } else {
    .check_d(x, y)
    X <- dplyr::bind_rows(
      "C1" = as.data.frame(x),
      "C2" = as.data.frame(y),
      .id = "type"
    )
    X$type <- as.factor(X$type)
  }
  min_size <- min(table(X$type))
  X <- X %>%
    dplyr::group_by(type) %>%
    dplyr::slice_sample(n = min_size) %>%
    dplyr::ungroup()
  training_set <- createDataPartition(X$type, p = split)
  ref <- caret::train(type ~ .,
                      data = X[training_set$Resample1, ],
                      method = method,
                      metric = "Accuracy",
                      trControl = control,
                      ...)
  test_res <- caret::predict.train(ref, newdata = X[-training_set$Resample1, ])
  p_hat <- mean(test_res == X[-training_set$Resample1, ]$type) - thresh
  min_accuracy <- max(table(X$type)) / nrow(X)
  pval <- stats::pbinom(p_hat * length(test_res),
                        size = length(test_res),
                        prob = min_accuracy,
                        lower.tail = FALSE)
  return(list("statistic" = p_hat, "p.value" = pval))
}
