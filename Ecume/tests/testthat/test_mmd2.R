library(testthat)

test_that("the mmd_test works with all inputs",{
  set.seed(08)
  x <- matrix(c(stats::runif(100, 0, 1),
                stats::runif(100, -1, 1)),
              ncol = 2)
  y <- matrix(c(stats::runif(100, 0, 1),
                stats::runif(100, -1, 1)),
              ncol = 2)
  set.seed(210)
  test <- mmd_test(x, y)
  expect_is(test, "list")
  set.seed(210)
  test2 <- mmd_test(y, x)
  expect_equal(test$statistic, test2$statistic)
  # linear
  set.seed(20)
  test <- mmd_test(x, y, type = "linear")
  expect_is(test, "list")
  expect_error(mmd_test(x, y, type = "linear", null = "exact"))
})
