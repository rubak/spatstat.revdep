## ----setup, include=FALSE, cache=FALSE---------------------------
library(knitr)
options(formatR.arrow=TRUE,width=67)

## ----------------------------------------------------------------
library(siplab)
trees <- as.data.frame(spruces)
summary(trees)

## ----------------------------------------------------------------
names(trees)[3] <- "dbh"
trees$dbh <- trees$dbh * 100
head(trees)

## ----------------------------------------------------------------
trees <- ppp(trees$x, trees$y, c(0,56), c(0,38), marks=trees$dbh)
summary(trees)

## ----------------------------------------------------------------
hegyi <- pairwise(trees, maxR=6, kernel=powers_ker,
    kerpar=list(pi=1, pj=1, pr=1, smark=1))

## ----------------------------------------------------------------
head(marks(hegyi))

## ----------------------------------------------------------------
hegyi_trim <- edges(hegyi, -6)
summary(hegyi_trim)

