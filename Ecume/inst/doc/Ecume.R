## ----packages, include=F------------------------------------------------------
library(knitr)
opts_chunk$set(
  fig.pos = "!h", out.extra = "", warning = FALSE, message = FALSE, 
  fig.align = "center", echo = FALSE
)
library(stats)

## -----------------------------------------------------------------------------
set.seed(20)
x <- rnorm(100, 0, 1)
y <- rnorm(200, 0, 1)
ks.test(x, y)

## -----------------------------------------------------------------------------
set.seed(20)
x <- rnorm(100, 0, 1)
y <- rnorm(200, 0, 2)
ks.test(x, y)

## -----------------------------------------------------------------------------
library(Ecume)
set.seed(20)
x <- rnorm(100, 0, 1)
w_x <- runif(100, 0, 1)
y <- rnorm(200, 0, 1)
w_y <- runif(200, 0, 1)
ks_test(x = x, y = y, w_x = w_x, w_y = w_y, thresh = .01)

## -----------------------------------------------------------------------------
set.seed(20)
x <- rnorm(100, 0, 1)
w_x <- runif(100, 0, 1)
y <- rnorm(200, 0, 2)
w_y <- runif(200, 0, 1)
ks_test(x = x, y = y, w_x = w_x, w_y = w_y, thresh = .01)

## -----------------------------------------------------------------------------
set.seed(20)
x <- matrix(c(rnorm(100, 0, 1),
              rnorm(100, 0, 1)),
            ncol = 2)
y <- matrix(c(rnorm(200, 0, 2),
              rnorm(200, 0, 1)),
            ncol = 2)
mmd_test(x = x, y = y, iterations = 10^4)

## -----------------------------------------------------------------------------
set.seed(20)
x <- matrix(c(rnorm(100, 0, 1),
              rnorm(100, 0, 1)),
            ncol = 2)
y <- matrix(c(rnorm(200, 0, 2),
              rnorm(200, 0, 1)),
            ncol = 2)
mmd_test(x = x, y = y, iterations = 10^4, type = "linear")

## -----------------------------------------------------------------------------
set.seed(20)
x <- matrix(c(rnorm(200, 0, 1),
              rnorm(200, 0, 1)),
            ncol = 2)
y <- matrix(c(rnorm(200, 0, 2),
              rnorm(200, 0, 1)),
            ncol = 2)
classifier_test(x = x, y = y)

## -----------------------------------------------------------------------------
set.seed(20)
x1 <- matrix(c(rnorm(200, 0, 1),
              rnorm(200, 0, 1)),
            ncol = 2)
x2 <- matrix(c(rnorm(200, 0, 2),
              rnorm(200, 0, 1)),
            ncol = 2)
x3 <- matrix(c(rnorm(200, 1, 1),
               rnorm(200, 0, 1)),
            ncol = 2)
classifier_test(x = list("x1" = x1, "x2" = x2, "x3" = x3))

