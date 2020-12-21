#
# util.R
#
#  $Revision: 1.17 $ $Date: 2017/02/04 05:59:09 $
#

gridproxy <- function(P, ..., dimyx=NULL, eps=NULL, xy=NULL, weights=NULL) {
  stopifnot(is.ppp(P))
  W <- as.owin(P)
  if(is.null(dimyx) && is.null(eps) && is.null(xy))
    dimyx <- 10
  M <- as.mask(W, dimyx=dimyx, eps=eps, xy=xy) 
  xy <- raster.xy(M, drop=TRUE)
  G <- as.ppp(xy, W=W)
  id <- nncross(G,P, what="which")
  if(!is.null(weights)) {
    # aggregate weights of P onto G
    check.nvector(weights, npoints(P))
    revid <- nncross(P, G, what="which")
    frevid <- factor(revid, levels=seq_len(npoints(G)))
    attr(id, "inverse") <- revid 
    attr(id, "weights") <- tapplysum(weights, list(frevid))
  }
  return(id)
}

FirstExtantEntry <- function(xlist, tags, whinge="No match") {
  y <- xlist[tags]
  y <- y[!unlist(lapply(y, is.null))]
  if(length(y) == 0)
    stop(whinge)
  return(names(y)[1])
}
  
sample.imagelist <- function(X, V) {
  Xvals <- lapply(X, safelookup, x=V)
  Xmat <- as.matrix(as.data.frame(Xvals))
  # ensure the entry names are not mangled by as.data.frame
  colnames(Xmat) <- names(X)
  return(Xmat)
}



applymaps <- local({
  
  applymaps <- function(maplist, x) {
    if(is.null(maplist)) return(x)
    if(is.language(maplist)) maplist <- list(maplist) else 
    stopifnot(is.list(maplist) && all(unlist(lapply(maplist, is.language))))
    if(is.data.frame(x)) {
      x <- as.data.frame(x)
      xenv <- list2env(as.list(x))
      y <- lapply(maplist, eval, envir=xenv)
      y <- as.data.frame(y)
    } else if(is.matrix(x)) {
      x <- as.data.frame(x)
      xenv <- list2env(as.list(x))
      y <- lapply(maplist, eval, envir=xenv)
      y <- as.matrix(as.data.frame(y))
    } else if(is.numeric(x)) {
      xenv <- list2env(as.list(x))
      y <- lapply(maplist, eval, envir=xenv)
      y <- unlist(y)
    } else if(is.list(x) && all(unlist(lapply(x, is.im))))  {
      xenv <- list2env(x)
      y <- lapply(maplist, evalim, xenv=xenv)
    } else if(inherits(x, "ssf")) {
      y <- ssf(unmark(x), applymaps(maplist, marks(x)))
    }
    return(y)
  }

  evalim <- function(f, xenv) eval(substitute(eval.im(ex, envir=en),
                                              list(ex=f, en=xenv)))

  applymaps
})
