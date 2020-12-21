## ----setup, include=FALSE, cache=FALSE---------------------------
library(knitr)
opts_chunk$set(fig.align='center', # fig.show='hold',
dev='pdf', out.width='.45\\textwidth') # , highlight=FALSE)
options(width=67)
library(siplab)

## ----------------------------------------------------------------
# distance from plant 1 at point (spot) x:
curve(abs(x - 2), from=0, to=10, ylim=c(0,4), lty=2,
    asp=1, ylab="distance")
curve(abs(x - 7), lty=2, add=T)  # distance from plant 2
curve(pmin(abs(x - 2), abs(x - 7)), add=T)  # minimum distance

## ----------------------------------------------------------------
# Assume siplab already installed (install.package("siplab")),
#   and loaded with library(siplab)
# The data must be in a marked ppp object:
threeTrees <- ppp(x=c(2,7,6), y=c(3,3,7), c(0,10), c(0,10),
    marks=c(10,10,10))  # marks are arbitrary (for now)
# Influence function. Takes distance components
cone_inf <- function(dx, dy, ...){  # and allow other args
    10 - sqrt(dx^2 + dy^2)  # 10 m height and radius at the base
}
# That's it, do it
a <- assimilation(threeTrees, influence=cone_inf)
points(a)  # add the tree locations to the influence map
# With a larger data set:
b <- assimilation(spruces, influence=cone_inf)

## ----------------------------------------------------------------
a
marks(a)
sum(marks(a)$aindex)

## ----------------------------------------------------------------
marks(threeTrees) <- c(35,30,40)  # size
# Cone, increased resolution for sharper boundaries:
a <- assimilation(threeTrees, pixsize=0.05)
# Paraboloid:
b <- assimilation(threeTrees, pixsize=0.05, infpar=c(a=2,
    b=.8, smark=1))

## ----------------------------------------------------------------
f <- function(x, size) {gnomon_inf(x, 0, size, par=c(a=2, b=1,
    smark=1))}
curve(f(x, 35), from=-1, to=6, lty=2, ylab="Influence") # tree 1
curve(f(5 - x, 30), lty=2, add=T)  # tree 2
curve(pmax(f(x, 35), f(5 - x, 30)), add=T)

## ----------------------------------------------------------------
curve(tass_inf(x, 0, marks=6), from=-3, to=3, asp=1, ylab="")
curve(gnomon_inf(x, 0, 6, par=c(a=1.3, b=2, smark=1)), lty=2,
    add=T)
# (for comparison)

## ----------------------------------------------------------------
a <- assimilation(finpines, afree=TRUE, influence=tass_inf,
    infpar=list(b=3.432, c=6.1, smark="height"))
aok <- edges(a, -2)  # remove trees near the plot border
head(marks(aok))

## ----------------------------------------------------------------
a <- assimilation(finpines, afree=TRUE, influence=tass_inf,
    infpar=list(b=3.432, c=6.1, smark="height"), asym=1)
head(marks(edges(a, -2)))

## ----------------------------------------------------------------
incr <- rnorm(30, 1, 0.2)  # diameter increments for 30 trees
incr
# Grow the diameter for 10 years:
D <- incr
D <- D + incr  # the long way, for clarity
D <- D + incr
D <- D + incr
D <- D + incr
D <- D + incr
D <- D + incr
D <- D + incr
D <- D + incr
D <- D + incr
regr <- lm(incr ~ D)  # regress increment over D at age 10
plot(incr ~ D)  # plot it
abline(regr)

## ----eval=FALSE--------------------------------------------------
#  trees <- boreasSA
#  dlim <- 3  # displacement limit
#  tolerance <- 0.1  # for convergence criterion
#  xy0 <- coords(trees)  # initial coordinates
#  lastdxy <- 0  # previous displacement
#  repeat {
#      a <- assimilation(trees, infpar=list(a=1, b=2.7,
#          smark="height"),centroid = TRUE)
#      dxy <- marks(a)[, c('cx','cy')] - xy0  # potential displ.
#      dxy[marks(a)$aindex <= 0,] <- 0  # ignore over-topped trees
#      d2 <- rowSums(dxy^2) # squared displacement lengths
#      toofar <- d2 > dlim^2
#      dxy[toofar, ] <- dlim * dxy[toofar, ] / sqrt(d2[toofar])
#      coords(trees) <- xy0 + dxy
#      if(max(abs(dxy - lastdxy)) < tolerance) break  # converged
#      lastdxy <- dxy
#  }
#  par(mfcol=1:2)
#  plot(edges(boreasSA, -5), main="Before", use.marks=F)
#  plot(edges(trees, -5), main="After", use.marks=F)

