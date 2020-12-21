# This one has leave-one-out option
expectedPairs_kernloo <- function(X, hx,hy, sigma=bw.CvL,
                    tol=.005, maxeval=1e6, maxsamp=5e3, leaveoneout=TRUE) {
    
    # Validate arguments
    stopifnot(inherits(X, "ppp"))
    # sort X by x coordinate
    ind <- order(X$x)
    coords(X) <- coords(X)[ind,]

    # take care of sigma
    if (is.function(sigma)) sigma <- sigma(X)
    stopifnot(is.numeric(sigma))

    stopifnot(is.numeric(tol) && tol > 0)

    # Check windows
    W <- as.owin(X)
    stopifnot(W$xrange==c(0,1) && W$yrange==c(0,1))

    cutoff <- cutoff2Dkernel("gaussian", sigma=sigma, varcov=NULL,
                           scalekernel=TRUE, cutoff=NULL, fatal=TRUE)

    # allocate results
    ep <- numeric(length(hx))
    ep2 <- numeric(length(hx))
    sampwts <- numeric(length(hx))

    # winwts account for the W \cap W_{-h} weighting of the window.
    # They are equal to \gamma if \rho_1 = \rho_2 = 1 everywhere
    if (is.rectangle(W)) {
        l <- diff(W$xrange)
        h <- diff(W$yrange)
        stopifnot(l == 1 && h == l)
        winwts <- (l - abs(hx))*(h - abs(hy))
    } else {
        #TODO: make them supported
        stop("non-rectangular windows not yet supported")
    }

    neghx <- hx < 0
    hx[neghx] <- -hx[neghx]
    hy[neghx] <- -hy[neghx]

    horder <- order(hx)
    allhx <- hx[horder]
    allhy <- hy[horder]

    #looping state
    valid_h <- which(winwts > 0)
    nloop <- 0
    inds <- valid_h
    nh_active <- length(valid_h)
    if (nh_active == 0) return(ep)

    se_from_sums <- function(s1, s2, n) sqrt( (n*s2/s1^2 - 1) / (n-1))

    # MC loop
    while (nh_active > 0) {
        weights <- sampwts[inds]
        hx <- allhx[inds]
        hy <- allhy[inds]

        minhx <- hx[1]

        # how many samples this iter?
        n_U <- max(1, min(floor(maxeval/nh_active), maxsamp))

        u <- sort(runif(n_U,0,1-minhx))
        v <- runif(n_U)

        f_at_u <- rhorho(u,v, hx,hy, X, sigma, cutoff, sorted=c("u","h","x"), leaveoneout=TRUE)

        weights <- weights + f_at_u$samps

        ep[inds] <- (ep[inds] + f_at_u$s)
        ep2[inds] <- (ep2[inds] + f_at_u$s2)

        sampwts[inds] <- weights
        sd_est <- (sqrt( (weights*ep2[inds] - ep[inds]^2) / (weights - 1)) / ep[inds])
        sd_est2 <- se_from_sums(ep[inds], ep2[inds], weights)

        passed <- sd_est < tol
        inds <- inds[!is.na(sd_est) & !passed]
        nh_active <- length(inds)
        nloop <- nloop + 1

        print(c(nh_active, min(weights), max(weights), min(sd_est), max(sd_est)), digits=2)
    }

    ep[valid_h] <- ep[valid_h] * winwts[valid_h] / sampwts[valid_h]
    ep[horder] <- ep
}

expectedPairs_iso_kernloo <- function(X, r, sigma=bw.CvL, tol=.001, maxeval=1e6,
                                    maxsamp=5e3, leaveoneout=TRUE) {

    # Validate arguments
    stopifnot(inherits(X, "ppp"))
    # sort X by x coordinate
    ind <- order(X$x)
    coords(X) <- coords(X)[ind,]

    if (is.function(sigma)) sigma <- sigma(X)
    stopifnot(is.numeric(sigma))

    stopifnot(is.numeric(tol) && tol > 0)

    # Check windows
    W <- as.owin(X)

    cutoff <- cutoff2Dkernel("gaussian", sigma=sigma, varcov=NULL,
                           scalekernel=TRUE, cutoff=NULL, fatal=TRUE)

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
        u <- U$x
        v <- U$y

        # sample 10 random directions
        for (j in 1:ndir) {
            # get a random direction
            hx <- rnorm(1)
            hy <- rnorm(1)
            rh <- sqrt(hx^2 + hy^2)
            hx <- hx/rh
            hy <- hy/rh

            f_at_u <- rhorho(u,v,hx*r,hy*r, X, sigma=sigma, cutoff=cutoff, leaveoneout=leaveoneout)

            weights <- weights + f_at_u$samps

            epr[inds] <- (epr[inds] + f_at_u$s)
            epr2[inds] <- (epr2[inds] + f_at_u$s2)
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
