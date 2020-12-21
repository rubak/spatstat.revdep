rrmix <- function(g, rhotspots, rsds, rweights, rbase = 1, log = TRUE) {
  if(!is.im(g)) stop("'g' must be of spatstat class 'im'")
  g <- g/integral(g)

  # w <- window
  # wg <- Window(g)
  # if(is.null(w)) w <- as.polygonal(wg)
  # if(!is.owin(w)) stop("'window' must be of spatstat class 'owin'")
  # w <- as.polygonal(w)
  # if(is.null(intersect.owin(wg,w,fatal=FALSE))) stop("'g' and 'window' must overlap")

  if(!is.matrix(rhotspots)||!is.numeric(rhotspots)) stop("'rhotspots' must be a numeric matrix")
  if(nrow(rhotspots)!=2) stop("'rhotspots' must have 2 rows")
  n <- ncol(rhotspots)

  if(!is.vector(rsds)||!is.numeric(rsds)||!is.vector(rweights)||!is.numeric(rweights)) stop("'rsds' and 'rweights' must be a numeric vector")
  if(length(rsds)!=n||length(rweights)!=n) stop("'rsds' and 'rweights' must have the same length as 'rhotspots' columns")
  if(any(rsds<=0)) stop("'rsds' must be strictly positive")

  rbase <- rbase[1]
  if(!is.numeric(rbase)) stop("'rbase' must be numeric")

  xy <- expand.grid(g$xcol,g$yrow)
  gdim <- dim(g)
  pgd <- prod(gdim)
  rvec <- rbase
  for(i in 1:n){
    hsi <- rhotspots[,i]
    ptx <- rep(hsi[1],pgd)
    pty <- rep(hsi[2],pgd)
    rvec <- rvec + rweights[i]*exp(-0.5/rsds[i]^2*((xy[,1]-ptx)^2+(xy[,2]-pty)^2))
  }
  rmat <- matrix(rvec,gdim[1],gdim[2],byrow=TRUE)
  gmat <- as.matrix(g)

  rint <- sum(rmat*gmat*g$xstep*g$ystep,na.rm=TRUE)
  rmat <- rmat/rint
  r <- im(rmat,xcol=g$xcol,yrow=g$yrow)[Window(g),drop=FALSE]
  f <- r * g
  if(log) r <- log(r)

  result <- solist(f=f,g=g,r=r)
  class(result) <- "rrim"
  return(result)
}
