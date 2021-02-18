library(testthat)

test_that("the wasserstein permutation test works",{
  set.seed(08)
  x <- matrix(c(stats::runif(20, 0, 1),
                stats::runif(20, -1, 1)),
              ncol = 2)
  y <- matrix(c(stats::runif(20, 0, 1),
                stats::runif(20, -1, 1)),
              ncol = 2)
  test <- wasserstein_permut(x, y)
  expect_is(test, "list")
  expect_true(test$statistic > 0)
  test_fast <- wasserstein_permut(x, y, fast = TRUE, S = 10)
  expect_is(test_fast, "list")
  expect_true(test_fast$statistic > 0)
  expect_error(wasserstein_permut(x, y, fast = TRUE))
  x <- matrix(c(stats::runif(20, 0, 1),
                stats::runif(20, -1, 1)),
              ncol = 1)
  y <- matrix(c(stats::runif(20, 0, 1),
                stats::runif(20, -1, 1)),
              ncol = 1)
  expect_is(wasserstein_permut(x, y), "list")
  expect_is(wasserstein_permut(x, y, fast = TRUE, S = 10), "list")
})
