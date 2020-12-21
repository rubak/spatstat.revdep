## Convert 'SpatialLines*' object to spatstat 'linnet' object
## 
## For 'SpatialLinesDataFrame', the data columns are copied 
## to the network as marks associated with the network segments.
##
## If fuse=TRUE, the code searches for pairs of points with the same (x,y)
##               coordinates that occur in different polylines,
##               and merges them together as identical vertices of the network.
##
##    Last edit: 2020/04/20 Adrian Baddeley 

if (!isClass("linnet"))
	setClass("linnet")

as.linnet.SpatialLines <- function(X, ..., fuse=TRUE) {
  if (!is.na(sp::is.projected(X)) && !sp::is.projected(X))
    stop("Only projected coordinates may be converted to spatstat class objects")
  if(!requireNamespace("spatstat", quietly = TRUE)) 
    stop("package spatstat is required for as.linnet.SpatialLines")
  #' extract bounding box to use as window
  bb <- bbox(X)
  BB <- spatstat::owin(bb[1,], bb[2,])
  #' 
  n <- length(X)
  xx <- yy <- numeric(0)
  ii <- jj <- integer(0)
  if(n > 0) {
    #' coordinates of all vertices
    crdlists <- coordinates(X)
    rowcounts <- unname(lapply(crdlists, function(x) sapply(x, nrow)))
    colcounts <- unname(lapply(crdlists, function(x) sapply(x, ncol)))
    if(any(unlist(colcounts) != 2)) stop("Coordinates should be 2-column matrices", call.=FALSE)
    #' 'rbind' all the matrices of coordinates
    xy <- unlist(lapply(crdlists, function(x) lapply(x, t)))
    xy <- matrix(xy, ncol=2, byrow=TRUE)
    xx <- xy[,1]
    yy <- xy[,2]
    #' construct indices for each level of list in original data
    ii <- rep(seq_along(crdlists), sapply(rowcounts, sum))
    jj <- unlist(lapply(rowcounts, function(x) rep(seq_along(x), as.integer(x))))
    #' check for *repeated* vertices within the same line
    rpt <- c(FALSE, (diff(xx) == 0) & (diff(yy) == 0) & (diff(ii) == 0) & (diff(jj) == 0))
    if(any(rpt)) {
      warning("Repeated vertices (on the same line) were removed", call.=FALSE)
      retain <- !rpt
      xx <- xx[retain]
      yy <- yy[retain]
      ii <- ii[retain]
      jj <- jj[retain]
    }
  }
  #' extract vertices 
  V <- spatstat::ppp(xx, yy, window=BB, check=!fuse)
  nV <- length(xx)
  #' join them
  edges <- NULL
  iii <- jjj <- integer(0)
  if(nV > 1) {
    seqn <- seq_len(nV)
    from <- seqn[-nV]
    to   <- seqn[-1]
    ok   <- diff(ii) == 0 & diff(jj) == 0
    from <- from[ok]
    to   <- to[ok]
    iii  <- ii[c(ok, FALSE)] #' indices backward
    jjj  <- jj[c(ok, FALSE)]
    if(fuse) {
      umap <- spatstat::uniquemap(V)
      retain <- (umap == seq_along(umap))
      V <- V[retain]
      renumber <- cumsum(retain)
      from <- renumber[umap[from]]
      to   <- renumber[umap[to]]
    }
    edges <- cbind(from, to)
  } 
  if(!is.null(edges)) {
    up <- (from < to)
    ee <- cbind(ifelse(up, from , to), ifelse(up, to, from))
    if(anyDuplicated(ee)) {
      u <- !duplicated(ee)
      from <- from[u]
      to   <- to[u]
      iii  <- iii[u]
      jjj  <- jjj[u]
    }
  }
  result <- spatstat::linnet(vertices=V, edges = edges, sparse=TRUE)
  if(spatstat::nsegments(result) == length(iii)) {
    df <- data.frame(LinesIndex=iii, LineIndex=jjj)
    if(.hasSlot(X, "data")) {
      DF <- slot(X, "data")
      df <- cbind(DF[iii,,drop=FALSE], df)
    }
    spatstat::marks(result$lines) <- df
  } else warning("Internal error: could not map data frame to lines")
  return(result)
}

setAs("SpatialLines", "linnet", function(from) as.linnet.SpatialLines(from))

setAs("SpatialLinesDataFrame", "linnet", function(from) as.linnet.SpatialLines(from))

