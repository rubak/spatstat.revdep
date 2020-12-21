## ------------------------------------------------------------------------
library("polyCub")

## ----example-f-----------------------------------------------------------
f <- function (s, sigma = 5)
{
    exp(-rowSums(s^2)/2/sigma^2) / (2*pi*sigma^2)
}

## ----example-polygon-----------------------------------------------------
hexagon <- list(
    list(x = c(7.33, 7.33, 3, -1.33, -1.33, 3),
         y = c(-0.5, 4.5, 7, 4.5, -0.5, -3))
)

## ----example, fig.width = 3, fig.height = 2.5----------------------------
plotpolyf(hexagon, f, xlim = c(-8,8), ylim = c(-8,8))

## ----product-Gauss, echo = -1, fig.show = "hold"-------------------------
par(mar = c(3,3,1,2))
polyCub.SV(hexagon, f, nGQ = 3, plot = TRUE)

## ------------------------------------------------------------------------
nrow(polyCub.SV(hexagon, f = NULL, nGQ = 3)[[1]]$nodes)

## ------------------------------------------------------------------------
polyCub.SVn <- function (polyregion, f, ..., nGQ = 20) {
    nw <- polyCub.SV(polyregion, f = NULL, ..., nGQ = nGQ)
    ## nw is a list with one element per polygon of 'polyregion'
    res <- sapply(nw, function (x)
        c(result = sum(x$weights * f(x$nodes, ...)), nEval = nrow(x$nodes)))
    structure(sum(res["result",]), nEval = sum(res["nEval",]))
}
polyCub.SVn(hexagon, f, nGQ = 3)

## ------------------------------------------------------------------------
for (nGQ in c(1:5, 10, 20)) {
    result <- polyCub.SVn(hexagon, f, nGQ = nGQ)
    cat(sprintf("nGQ = %2i: %.12f (n=%i)\n", nGQ, result, attr(result, "nEval")))
}

## ---- message = FALSE----------------------------------------------------
library("spatstat")
hexagon.owin <- owin(poly = hexagon)

## ----midpoint, echo = -1, fig.show = "hold"------------------------------
par(mar = c(3,3,1,3), xaxs = "i", yaxs = "i")
polyCub.midpoint(hexagon.owin, f, eps = 0.5, plot = TRUE)

## ------------------------------------------------------------------------
intrfr <- function (R, sigma = 5)
{
    (1 - exp(-R^2/2/sigma^2))/2/pi
}

## ------------------------------------------------------------------------
polyCub.iso(hexagon, intrfr = intrfr, center = c(0,0))

## ------------------------------------------------------------------------
gpclibPermit()  # accept gpclib license (prohibits commercial use)
polyCub.exact.Gauss(hexagon.owin, mean = c(0,0), Sigma = 5^2*diag(2))

