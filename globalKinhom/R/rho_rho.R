rhorho <- function(ux, uy, hx, hy, X, sigma, cutoff=8*sigma, sorted=c(),leaveoneout=TRUE) {

    if (!any("x" == sorted)) {
        inds <- order(X$x)
        coords(X) <- coords(X)[inds,]
    }
    x <- X$x
    y <- X$y

    if (!any("u" == sorted)) {
        indsq <- order(ux)
    } else {
        indsq <- 1:length(ux)
    }
    xq <- ux[indsq]
    yq <- uy[indsq]

    if (!any("h" == sorted)) {
        indsh <- order(hx)
        hx <- hx[indsh]
        hy <- hy[indsh]
    } else {
        indsh <- 1:length(hx)
    }

    wq <- diggle_weights(ux, uy, sigma)

    up <- outer(ux, hx, `+`)
    vp <- outer(uy, hy, `+`)
    inside <- inside.owin(up, vp, as.owin(X))
    wp <- matrix(0,nrow=length(ux), ncol=length(hx))
    wp[inside] <- diggle_weights(up[inside], vp[inside], sigma)

    if (leaveoneout) {
        zz <- .C("rho_rho_excess",
            nquery = as.integer(length(ux)),
            xq = as.double(xq),
            yq = as.double(yq),
            ndata = as.integer(length(x)),
            xd = as.double(x),
            yd = as.double(y),
            nsep = as.integer(length(hx)),
            xh = as.double(hx),
            yh = as.double(hy),
            rmaxi = as.double(cutoff),
            sig= as.double(sigma),
            result = as.double(double(length(ux)*length(hx))))
            #package = "globalKinhom")

        res <- zz$result
        dim(res) <- c(length(hx), length(ux))
        res <- t(res) * wq * wp
    } else res <- matrix(0,nrow=length(ux),ncol=length(hx))

    dfun <- densityfun(X,sigma)
    rhou <- dfun(ux,uy) 

    rhoup <- matrix(0,nrow=length(ux), ncol=length(hx))
    rhoup[inside] <- dfun(up[inside], vp[inside])

    res[, indsh] <- (rhou * rhoup) - res
    inside[,indsh] <- inside

    data.frame(s=colSums(res), s2=colSums(res^2), samps=colSums(inside))
    
}

crossrhorho <- function(ux,uy,hx,hy,X,Y,sigma) {

    dfun <- densityfun(X,sigma)
    dfun2 <- densityfun(Y,sigma)
    rhou <- dfun(ux,uy) 

    up <- outer(ux, hx, `+`)
    vp <- outer(uy, hy, `+`)

    inside <- inside.owin(up,vp,as.owin(X))

    rhoup <- matrix(0,nrow=length(ux), ncol=length(hx))
    rhoup[inside] <- dfun2(up[inside], vp[inside])

    data.frame(s=colSums(rhou * rhoup), s2=colSums(rhou^2 * rhoup^2), samps=colSums(inside))
}
