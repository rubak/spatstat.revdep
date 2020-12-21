###############################################################################
# These are numerical estimates of the quatities f and f_12 in the paper.
# Perhaps a better name would be "expected number of pairs for poisson process(es)".
# Hence, naming these functions expectedPairs and the like

# expectedPairs is just a wrapper for expectedCrossPairs that takes the right args
expectedPairs <- function(rho, hx, hy=NULL, method=c("mc", "lattice"), tol=.005,
                            dx=diff(as.owin(rho)$xrange)/200, maxeval=1e6,
                            maxsamp=5e3) {
    expectedCrossPairs(rho1=rho, rho2=NULL, hx, hy, method, tol, dx,
                        maxsamp=maxsamp, maxeval=maxeval)
}

expectedCrossPairs <- function(rho1, rho2=NULL, hx, hy=NULL, method=c("mc", "lattice"),
                            tol=.005, dx=diff(as.owin(rho1)$xrange)/200,
                            maxeval=1e6, maxsamp=5e3) {

    # Validate arguments
    stopifnot(is.im(rho1) || inherits(rho1, "funxy"))
    if (is.im(rho1)) {
        rho1.im <- rho1
        rho1 <- funxy(function(x,y) interp.im(rho1.im, x, y), as.owin(rho1.im))
    }
    cross <- !is.null(rho2)
    if (cross) {
        stopifnot(is.im(rho2) || inherits(rho2, "funxy"))
        if (is.im(rho2)) {
            rho2.im <- rho2
            rho2 <- funxy(function(x,y) interp.im(rho2.im, x, y), as.owin(rho2.im))
        }
    } else {
        rho2 <- rho1 # for convenience later
    }
    method <- match.arg(method)

    use.mc <- FALSE
    use.lattice <- FALSE
    if (method == "mc") {
        use.mc <- TRUE
        stopifnot(is.numeric(tol) && tol > 0)
    } else if (method == "lattice") {
        use.lattice <- TRUE
        stopifnot(is.numeric(dx) && dx > 0 && dx < as.owin(rho1))
    }

    # Check windows...
    W <- as.owin(rho1)
    if (cross) {
        W2 <- as.owin(rho2)
        stopifnot(all(W$xrange == W2$xrange) && all(W$yrange == W2$yrange))
    }

    # Get coordinates of h
    xy <- xy.coords(hx, hy)
    allhx <- xy$x
    allhy <- xy$y
    nh <- length(hx)

    # allocate results
    epr <- numeric(length(hx))
    epr2 <- numeric(length(hx))
    sampwts <- numeric(length(hx))

    if (is.rectangle(W)) {
        winwts <- (abs(diff(W$xrange)) - abs(allhx))*(abs(diff(W$yrange)) - abs(allhy))
    } else {
        winwts <- numeric(length(allhx))
        for (i in 1:length(allhx)) {
            W2 <- shift(W, - c(allhx[i], allhy[i]))

            winwts[i] <- overlap.owin(W, W2)
        }
    }

    # state of looping
    valid_h <- which(winwts > 0) # if winwts<=0 then f(h) = 0
    nloop <- 0
    inds <- valid_h # hs that haven't converged
    nh_active <- length(valid_h)
    if (nh_active == 0) return(epr) # none of the hs are valid

    # Main loop
    while (nh_active > 0) {
        weights <- sampwts[inds]
        hx <- allhx[inds]
        hy <- allhy[inds]

        # how many samples to take this iteration?
        n_samp <- max(1, min(floor(maxeval/nh_active), maxsamp))

        # Monte carlo sample of points in W
        U <- runifpoint(n_samp, W)

        # Corresponding intensity(s)
        rho1U <- rho1(U)
        rho2U <- if (cross) rho2(U) else rho1U

        # Locations of second point, for each h
        Uplushx <- outer(U$x, hx, `+`)
        Uplushy <- outer(U$y, hy, `+`)
        Uminushx <- outer(U$x, hx, `-`)
        Uminushy <- outer(U$y, hy, `-`)

        # Are these in the window?
        Uplus_inside <- inside.owin(Uplushx, Uplushy, W)
        Uminus_inside <- inside.owin(Uminushx, Uminushy, W)

        # Where yes, corresponding intensity(s)
        rhoUplus <- matrix(0, nrow=n_samp, ncol=nh_active)
        rhoUminus <- matrix(0, nrow=n_samp, ncol=nh_active)
        rhoUplus[Uplus_inside] <- rho2(Uplushx[Uplus_inside], Uplushy[Uplus_inside])
        rhoUminus[Uminus_inside] <- rho1(Uminushx[Uminus_inside], Uminushy[Uminus_inside])

        # Number of new points for each h?
        weights <- weights + colSums(Uplus_inside + Uminus_inside)
        sampwts[inds] <- weights

        # Update estimates
        epr[inds] <- (epr[inds]
                        + (rho1U %*% rhoUplus)
                        + (rho2U %*% rhoUminus))

        epr2[inds] <- (epr2[inds]
                        + ((rho1U^2) %*% (rhoUplus^2))
                        + ((rho2U^2) %*% (rhoUminus^2)))

        # estimate the standard error of the estimates
        sd_est <- ( sqrt( (weights*epr2[inds] - epr[inds]^2) / (weights - 1))
                            / epr[inds])

        passed <- sd_est < tol
        inds <- inds[!passed]
        nh_active <- length(inds)
        nloop <- nloop + 1


        print(c(nh_active, min(weights), max(weights),
                    min(sd_est), max(sd_est)), digits=2)
    }
    # Final output. sampwts is the number of applicable monte carlo samples
    # winwts is the area of W \cap W_{-h}
    epr[valid_h] <- epr[valid_h] * winwts[valid_h] / sampwts[valid_h]
    epr
}

expectedPairs_iso <- function(rho, r, tol=.001, maxeval=1e6, maxsamp=5e3) {
    expectedCrossPairs_iso(rho, NULL, r, tol, maxeval, maxsamp)
}

expectedCrossPairs_iso <- function(rho1, rho2=NULL, r,
                    tol=.001, maxeval=1e6, maxsamp=5e3) {

    # Validate arguments
    stopifnot(is.im(rho1) || inherits(rho1, "funxy"))
    if (is.im(rho1)) {
        rho1.im <- rho1
        rho1 <- funxy(function(x,y) interp.im(rho1.im, x,y), as.owin(rho1.im))
    }
    cross <- !is.null(rho2)
    if (cross) {
        stopifnot(is.im(rho2) || inherits(rho2, "funxy"))
        if (is.im(rho2)) {
            rho2.im <- rho2
            rho2 <- funxy(function(x,y) interp.im(rho2.im, x,y), as.owin(rho2.im))
        }
    } else {
        rho2 <- rho1
    }

    stopifnot(is.numeric(tol) && tol > 0)

    # Check windows
    W <- as.owin(rho1)
    if (cross) {
        W2 <- as.owin(rho2)
        stopifnot(all(W$xrange == W2$xrange) && all(W$yrange== W2$yrange))
    }

    # allocate results
    epr <- numeric(length(r))
    epr2 <- numeric(length(r))
    sampwts <- numeric(length(r))

    # winwts account for the W \cap W_{-h} weighting of the window.
    # They are equal to \gamma_r if \rho_1 = \rho_2 = 1 everywhere
    if (is.rectangle(W)) {
        l <- diff(W$xrange)
        h <- diff(W$yrange)
        winwts <- (2*pi - 4*(r/l + r/h) + 2*r^2/(l*h)) * l*h/(2*pi)
    } else {
        #TODO: make them supported
        stop("non-rectangular windows not yet supported")
    }

    #looping state
    allr <- r
    valid_r <- which(winwts > 0)
    nloop <- 0
    inds <- valid_r
    nr_active <- length(valid_r)
    if (nr_active == 0) return(epr)

    se_from_sums <- function(s1, s2, n) sqrt( (n*s2/s1^2 - 1) / (n-1))

    # MC loop
    ndir <- 1
    while (nr_active > 0) {
        weights <- sampwts[inds]
        r <- allr[inds]

        # how many samples this iter?
        n_U <- max(1, min(floor(maxeval/nr_active), maxsamp))

        U <- runifpoint(n_U, W)
        rho1U <- rho1(U)

        # sample 10 random directions
        for (j in 1:ndir) {
            # get a random direction
            hx <- rnorm(1)
            hy <- rnorm(1)
            rh <- sqrt(hx^2 + hy^2)
            hx <- hx/rh
            hy <- hy/rh

            Uphx <- outer(U$x, hx*r, `+`)
            Uphy <- outer(U$y, hy*r, `+`)

            inside <- inside.owin(Uphx, Uphy, W)

            rhoUp <- matrix(0, nrow=n_U, ncol=nr_active)
            rhoUp[inside] <- rho2(Uphx[inside], Uphy[inside])

            weights <- weights + colSums(inside)

            epr[inds] <- (epr[inds] + (rho1U %*% rhoUp))
            epr2[inds] <- (epr2[inds] + (rho1U)^2 %*% (rhoUp)^2)
        }

        sampwts[inds] <- weights
        sd_est <- sqrt(ndir)*(sqrt( (weights*epr2[inds] - epr[inds]^2) / (weights - 1)) / epr[inds])

        passed <- sd_est < tol
        inds <- inds[!is.na(sd_est) & !passed]
        nr_active <- length(inds)
        nloop <- nloop + 1

        #print(c(nr_active, min(weights), max(weights), min(sd_est), max(sd_est)), digits=2)
    }

    epr[valid_r] <- epr[valid_r] * winwts[valid_r] / sampwts[valid_r]
    epr
}

diggle_weights <- function(x, y, sigma) { #TODO: support windows that aren't unit square?
#    x <- coords(X)$x
#    y <- coords(X)$y
#    W <- as.owin(X)
#    stopifnot("rectangle" == W$type)
#    xrange <- W$xrange
#    yrange <- W$yrange
    xrange <- c(0,1)
    yrange <- c(0,1)
    1/((pnorm((xrange[2]-x)/sigma) - pnorm((xrange[1] - x)/sigma)) *
                (pnorm((yrange[2] - y)/sigma) - pnorm((yrange[1] - y)/sigma)))
}
