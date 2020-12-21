edge.Iso <- function(X,Y=X,win=Window(X), paired=FALSE, dx=NULL, dy=NULL, d=NULL) {
    if (is.null(d)) {
        if (is.null(dx) && is.null(dy)) {
            if (paired) {
                dx = Y$x - X$x
                dy = Y$y - X$y
            } else {
                dx = outer(Y$x, X$x, `-`)
                dy = outer(Y$y, X$y, `-`)
            }
        }

        d <- sqrt(dx^2 + dy^2)
    }

        
    w <- switch(win$type,
        rectangle= {
                    width <- diff(win$xrange)
                    height <- diff(win$yrange)
                    (2*pi - 4*d*(1/width + 1/height) + 2*d^2/(width*height)) * width*height/(2*pi)
        })

    if (is.null(w)) stop("edge.Iso: unknown window type")

    w
}
