###############################################################################
# K_global from the paper
# argument list is based on spatstat Kinhom, though many of those arguments
# are not relevant here. ... are passed to density.ppp or densityfun.ppp
# if no lambda, use a kernel version of expectedCrossPairs
Kglobal <-
function(X, lambda=NULL, ..., sigma=bw.CvL(X), r=NULL, rmax=NULL, breaks=NULL,
            normtol=.005, discrete.lambda=FALSE,
            interpolate=TRUE, interpolate.fac=10, isotropic=TRUE,
            leaveoneout=TRUE, exp_prs=NULL,
            interpolate.maxdx=diameter(as.owin(X))/100, dump=FALSE) {
    # Check inputs
    verifyclass(X, "ppp")
    W <- as.owin(X)
    npts <- npoints(X)
    AreaW <- area(W)

    # What should the rs be?
    rfixed <- !missing(r) || !missing(breaks)
    rmaxdefault <- if (!is.null(rmax)) rmax else rmax.rule("K", W, npts/AreaW)
    breaks <- handle.r.b.args(r, breaks, W, rmaxdefault=rmaxdefault)
    r <- breaks$r
    rmax <- breaks$max

    lambda.given <- !is.null(lambda)
    ep.given <- !is.null(exp_prs)

    pairs <- closepairs(X, rmax, twice=FALSE)
    hx <- pairs$dx #pairs$xi - pairs$xj
    hy <- pairs$dy #pairs$yi - pairs$yj
    pairdist <- pairs$d

    # Get gammas, depending on isotropic and interpolate options
    if (isotropic) {
        if (interpolate) {
            dr <- min(sigma/interpolate.fac, interpolate.maxdx)
            rcheck <- seq(0, max(pairdist) + dr, by=dr)
        } else {
            rcheck <- pairdist
        }

        if (ep.given) {
            gamma <- exp_prs(rcheck)
        } else if (lambda.given) {
            gamma <- expectedPairs_iso(lambda,rcheck, tol=normtol)
        } else {
            gamma <- expectedPairs_iso_kernloo(X,rcheck,sigma=sigma,tol=normtol,leaveoneout=leaveoneout)
        }

        if (interpolate) {
            gamma <- approx(rcheck, gamma, pairdist, rule=2)$y
        }
    } else { # !isotropic
        if (interpolate) {
            dhx <- min(sigma/interpolate.fac, interpolate.maxdx)
            npt <- ceiling(rmax/dhx)
            xs <- (-npt:npt)*dhx
            lathx <- outer(xs, xs, function(x,y) x)
            lathy <- outer(xs, xs, function(x,y) y)
        } else {
            lathx <- hx
            lathy <- hy
        }

        if (ep.given) {
            gamma <- exp_prs(lathx, lathy)
        } else if (lambda.given) {
            gamma <- expectedPairs_iso(lambda,lathx,lathy,tol=normtol)
        } else {
            gamma <- expectedPairs_kernloo(X,lathx,lathy,sigma=sigma,tol=normtol,leaveoneout=leaveoneout)
        }

        if (interpolate) {
            dim(gamma) <- c(2*npt + 1, 2*npt + 1)
            latgamma.im <- as.im(t(gamma), xrow=xs, ycol=xs)

            gamma <- interp.im(latgamma.im, hx, hy)
        }
    }

    bins <- .bincode(pairdist, r, include.lowest=TRUE)
    K <- numeric(length(r))

    for (i in 2:length(r)) {
        K[i] <- sum(1/gamma[bins == i - 1])
    }
    K <- cumsum(K)*2 # make up for twice=FALSE above
    Kf <- data.frame(r=r, theo=pi*r^2, global=K)

    out <- fv(Kf, argu="r", ylab=quote(K(r)), valu="global", fmla= . ~ r,
        alim=c(0,rmax), labl=c("r", "%s[Pois](r)", "%s[global](r)"),
        desc=c("distance argument r", "theoretical poison %s",
        "global correction %s"), fname="K")

    if (dump) {
        attr(out, "gammas") <- gamma
        attr(out, "prs") <- pairs
    }

    out
}

# If lambdaX and lambdaY are both NULL, use the kernel version
# otherwise, use Monte Carlo estimates for f
Kcross.global <-
function(X, Y, lambdaX=NULL, lambdaY=NULL, ..., sigma=bw.CvL(X), r=NULL,
            rmax=NULL, breaks=NULL, normtol=.005,
            discrete.lambda=FALSE, interpolate=TRUE, isotropic=TRUE,
            interpolate.fac=10, leaveoneout=TRUE, exp_prs=NULL,
            interpolate.maxdx=diameter(as.owin(X))/100, dump=FALSE) {
    # Check inputs
    verifyclass(X, "ppp")
    verifyclass(Y, "ppp")
    W <- as.owin(X)
    W2 <- as.owin(Y)
    stopifnot(all(W$xrange == W2$xrange) && all(W$yrange == W2$yrange))

    npts <- npoints(X)
    AreaW <- area(W)

    # What should the rs be?
    rfixed <- !missing(r) || !missing(breaks)
    rmaxdefault <- if (!is.null(rmax)) rmax else rmax.rule("K", W, npts/AreaW)
    breaks <- handle.r.b.args(r, breaks, W, rmaxdefault=rmaxdefault)
    r <- breaks$r
    rmax <- breaks$max

    # How to compute gamma?
    ep.given <- !is.null(exp_prs)

    if (!ep.given) {
        if (!is.null(sig_tmp <- attr(lambdaX, "sigma"))) {
            sigma <- sig_tmp
        }

        lambdaX <- fixLambda(lambdaX, X, discrete.lambda, sigma, ...)
        lambdaY <- fixLambda(lambdaY, Y, discrete.lambda, sigma, ...)
    }

    pairs <- crosspairs(X, Y, rmax)
    hx <- pairs$dx #pairs$xi - pairs$xj
    hy <- pairs$dy #pairs$yi - pairs$yj
    pairdist <- pairs$d

    # Get gammas, depending on isotropic and interpolate options
    if (isotropic) {
        rh <- sqrt(hx^2 + hy^2)
        if (interpolate) {
            dr <- min(sigma/interpolate.fac, interpolate.maxdx)
            rcheck <- seq(0, max(rh) + dr, by=dr)
        } else {
            rcheck <- rh
        }

        if(ep.given) {
            if (is.function(exp_prs)) {
                gamma <- exp_prs(rcheck)
            } else {
                stop("exp_prs is unknown format")
            }
        } else {
            gamma <- expectedCrossPairs_iso(lambdaX, lambdaY, rcheck, tol=normtol)
        }

        if (interpolate) {
            gamma <- approx(rcheck, gamma, pairdist, rule=2)$y
        }
    } else { # !isotropic
        if (interpolate) {
            dhx <- min(sigma/interpolate.fac, interpolate.maxdx)
            npt <- ceiling(rmax/dhx)
            xs <- (-npt:npt)*dhx
            lathx <- outer(xs, xs, function(x,y) x)
            lathy <- outer(xs, xs, function(x,y) y)
        } else {
            lathx <- hx
            lathy <- hy
        }

        if (ep.given) {
            if (is.function(exp_prs)) {
                gamma <- exp_prs(lathx, lathy)
            } else {
                stop("exp_prs is unknown format")
            }
        } else {
            gamma <- expectedCrossPairs(lambdaX, lambdaY, lathx, lathy, tol=normtol)
        }

        if (interpolate) {
            dim(gamma) <- c(2*npt + 1, 2*npt + 1)
            latgamma.im <- as.im(t(gamma), xcol=xs, yrow=xs)

            gamma <- interp.im(latgamma.im, hx, hy)
        }
    }

    bins <- .bincode(pairdist, r, include.lowest=TRUE)
    K <- numeric(length(r))

    for (i in 2:length(r)) {
        K[i] <- sum(1/gamma[bins == i - 1])
    }

    K <- cumsum(K)

    Kf <- data.frame(r=r, theo=pi*r^2, global=K)

    out <- fv(Kf, argu="r", ylab=quote(K(r)), valu="global", fmla= . ~ r,
        alim=c(0,rmax), labl=c("r", "%s[Pois](r)", "%s[global](r)"),
        desc=c("distance argument r", "theoretical poison %s",
        "global correction %s"), fname="K")

    if (dump) {
        attr(out, "gammas") <- gamma
        attr(out, "prs") <- pairs
    }
    out
}


fixLambda <- function(lambdaX, X, discrete.lambda, sigma, ...) {
    W <- as.owin(X)
    if (is.null(lambdaX)) {
        if (discrete.lambda) {
            lambdaX.im <- density(X, sigma, ...)
            lambdaX <- funxy(function(x,y) interp.im(lambdaX.im, x,y), W)
        } else {
            lambdaX <- densityfun(X, sigma, ...)
        }
    } else {
        Wl <- as.owin(lambdaX)
        stopifnot(all(Wl$xrange == W$xrange) && all(Wl$yrange == W$yrange))
    }

    lambdaX
}
