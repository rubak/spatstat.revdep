library(testthat)

test_that("stouffer_Z_score works with all inputs",{
  set.seed(90)
  pvals <- stats::runif(100, 0, 1)
  weights <- stats::runif(100, 0, 1)
  expect_is(stouffer_zscore(pvals, weights), "list")
  expect_equal(stouffer_zscore(pvals, weights, side = 'right'),
               stouffer_zscore(1 - pvals, weights, side = 'left'))
  expect_equal(stouffer_zscore(1 - pvals, weights, side = 'right'),
               stouffer_zscore(pvals, weights, side = 'left'))
  pval <- stats::runif(1, 0, 1)
  combi <- stouffer_zscore(pval, side = 'right')
  expect_equal(pval, combi$p.value)
  pvals <- rep(pval, 100)
  for (i in 2:100) {
    combi <- stouffer_zscore(pvals[1:i], weights[1:i])
    expect_true(pval > combi$p.value)
  }
})
