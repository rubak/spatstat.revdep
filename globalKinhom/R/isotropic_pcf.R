# isotropic estimate for pcf. bw is the bandwidth for pcf estimation, _NOT_ for
# intensity estimation. ... are passed to densityfun.ppp()
pcfglobal <-
function(X, lambda=NULL, ..., sigma=bw.CvL(X), r=NULL, rmax=NULL,
    kernel="epanechnikov", bw=NULL, stoyan=0.15, normtol=.005, ratio=FALSE,
    discrete.lambda=FALSE, divisor=c("r", "d"),
    leaveoneout=TRUE, interpolate=TRUE, interpolate.fac=10, exp_prs=NULL,
    interpolate.maxdx=diameter(as.owin(X))/100, dump=FALSE) {

    verifyclass(X, "ppp")
    W <- as.owin(X)
    areaW <- area(W)
    npts <- npoints(X)

    rmaxdefault <- if (is.null(rmax)) rmax <- rmax.rule("K", W, npts/areaW)
    breaks <- handle.r.b.args(r, NULL, W, rmaxdefault=rmax)
    r <- breaks$r
    rmax <- breaks$max

    alim <- c(0, min(rmax, rmaxdefault))

    kernel <- match.kernel(kernel)

    if (is.null(bw) && (kernel == "epanechnikov")) {
        h <- stoyan/sqrt(npts/areaW)
        hmax <- h
        bw <- h/sqrt(5)
    }
    else if (is.numeric(bw)) {
        hmax <- 3 * bw
    }
    else {
        hmax <- 2 * stoyan/sqrt(npts/areaW)
    }
    denargs <- list(kernel=kernel, bw=bw, n=length(r), from=0, to=rmax)

    lambda.given <- !is.null(lambda)
    ep.given <- !is.null(exp_prs)

    if (interpolate) {
        dr <- min(sigma/ interpolate.fac, interpolate.maxdx)
        if (dr < r[2]) interpolate=FALSE
    }
    if (interpolate) {
        r_test <- seq(0, rmax - r[2] + dr, by=dr)
    } else {
        r_test <- r
    }

    if (ep.given) {
        gammas <- exp_prs(r_test)
    } else if (lambda.given) {
        gammas <- expectedPairs_iso(lambda, r_test, tol=normtol)
    } else {
        gammas <- expectedPairs_iso_kernloo(X, r_test, sigma=sigma, tol=normtol, leaveoneout=leaveoneout)
    }

    if (interpolate) {
        gammas <- approx(r_test, gammas, r, rule=2)$y
#         spl <- smooth.spline(r_test, gammas, df=length(r_test))
#         gammas <- predict(spl, r)$y
    }

    prs <- closepairs(X, rmax + hmax, what='ijd')

    df <- data.frame(r=r, theo=rep.int(1, length(r)))
    out <- ratfv(df, NULL, gammas, "r", quote(g(r)), "theo",
                NULL, alim, c("r", "%s[Pois](r)"), c("distance argument r",
                "theoretical Poisson %s"), fname="g", ratio=ratio)

    bw.used <- NULL

    kdenG <- sewpcf(prs$d, 1, denargs, 1, divisor=divisor)
    gG <- kdenG$g/gammas
    bw.used <- attr(kdenG, "bw")
    if (!ratio) {
        out <- bind.fv(out, data.frame(global=gG), "hat(%s)[global](r)",
                        "Global intensity reweighted estimate of %s", "global")
    } else {
        out <- bind.ratfv(out, data.frame(global=gG * gammas), gammas,
                        "hat(%s)[global](r)",
                        "Global intensity reweighted estimate of %s", "global")
    }

    formula(out) <- . ~ r
    fvnames(out, ".") <- setdiff(rev(colnames(out)), c("r", "v"))
    unitname(out) <- unitname(X)
    if (ratio)
        out <- conform.ratfv(out)
    attr(out, "bw") <- bw.used

    if (dump) {
        attr(out, "prs") <- prs
        attr(out, "fs") <- gammas
    }

    out
}

pcfcross.global <- function(X,Y, lambdaX=NULL, lambdaY=NULL, ...,
    sigma=bw.CvL(X), r=NULL, rmax=NULL, kernel="epanechnikov", bw=NULL,
    stoyan=0.15, normtol=.005, ratio=FALSE, discrete.lambda=FALSE,
    divisor=c("r", "d"), analytical=NULL, interpolate=TRUE,
    interpolate.fac=10, exp_prs=NULL,
    interpolate.maxdx=diameter(as.owin(X))/100, dump=FALSE) {

    verifyclass(X, "ppp")
    verifyclass(Y, "ppp")
    W <- as.owin(X)
    W2 <- as.owin(Y)
    stopifnot(all(W$xrange == W2$xrange) && all(W$yrange == W2$yrange))
    areaW <- area(W)
    npts <- npoints(X)

    rmaxdefault <- if (is.null(rmax)) rmax <- rmax.rule("K", W, npts/areaW)
    breaks <- handle.r.b.args(r, NULL, W, rmaxdefault=rmax)
    r <- breaks$r
    rmax <- breaks$max

    alim <- c(0, min(rmax, rmaxdefault))

    kernel <- match.kernel(kernel)

    # Bandwidth for the pcf estimate
    if (is.null(bw) && (kernel == "epanechnikov")) {
        h <- stoyan/sqrt(npts/areaW)
        hmax <- h
        bw <- h/sqrt(5)
    }
    else if (is.numeric(bw)) {
        hmax <- 3 * bw
    }
    else {
        hmax <- 2 * stoyan/sqrt(npts/areaW)
    }
    denargs <- list(kernel=kernel, bw=bw, n=length(r), from=0, to=rmax)

    ep.given <- !is.null(exp_prs)

    if (interpolate) {
        dr <- min(sigma/ interpolate.fac, interpolate.maxdx)
        if (dr < r[2])
            interpolate=FALSE
    }
    if (interpolate) {
        r_test <- seq(0, rmax - r[2] + dr, by=dr)
    } else {
        r_test <- r
    }

    if (ep.given){
        if (is.function(exp_prs)) {
            gammas <- exp_prs(r_test)
        } else {
            stop("exp_prs given in unknown form")
        }
    } else {
        lX <- fixLambda(lambdaX, X, discrete.lambda, sigma, ...)
        lY <- fixLambda(lambdaY, Y, discrete.lambda, sigma, ...)

        gammas <- expectedCrossPairs_iso(lX, lY, r_test, tol=normtol)
    }

    if (interpolate) {
        gammas <- approx(r_test, gammas, r, rule=2)$y
    }

    prs <- crosspairs(X, Y, rmax + hmax, what='all')

    df <- data.frame(r=r, theo=rep(1, length(r)))
    out <- ratfv(df, NULL, gammas, "r", quote(c(r)), "theo", NULL, alim,
                c("r", "%s[Pois](r)"),
                c("distance argument r", "theoretical poisson %s"),
                fname="c", ratio=ratio)

    bw.used <- NULL

    kdenG <- sewpcf(prs$d, 1, denargs, 1, divisor=divisor)
    cG <- kdenG$g/gammas
    bw.used <- attr(kdenG, "bw")
    if (!ratio) {
        out <- bind.fv(out, data.frame(global=cG), "hat(%s)[global](r)",
                "Global intensity reweighted estimate of %s", "global")
    } else {
        out <- bind.ratfv(out, data.frame(global=cG*gammas), gammas,
                "hat(%s)[global](r)", "Global intensity reweighted estimate of %s",
                "global")
    }

    formula(out) <- . ~ r
    fvnames(out, ".") <- setdiff(rev(colnames(out)), c("r", "v"))
    unitname(out) <- unitname(X)
    if (ratio)
        out <- conform.ratfv(out)

    attr(out, "bw") <- bw.used
    if (dump) {
        attr(out, "prs") <- prs
        attr(out, "fs") <- gammas
    }

    out
}
