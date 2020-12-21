sgmix <- function(mean, vcv, window, p0=0, p=NULL, resolution=128, int=1) {
  if(!is.matrix(mean)) stop("'mean' must be a matrix")
  if(nrow(mean)!=2) stop("'mean' is of incorrect dimension")
  n <- ncol(mean)

  if(int<=0) stop("'int' must be positive")
  if(resolution<=1) stop("'resolution' must be >= 1")

  if(!is.owin(window)) stop("'window' must be of spatstat class 'owin'")

  if(length(p0)>1) p0 <- p0[1]
  if(!is.numeric(p0)) stop("'p0' must be numeric")
  if(p0<0||p0>1) stop("'p0' must be in [0,1]")

  if(is.null(p)) p <- rep((1-p0)/n,n)
  if(!is.numeric(p)) stop("'p' must be numeric")
  if(length(p)!=n)
  if(any(p<0)||any(p>1)) stop("all elements of 'p' must be in [0,1]")
  if(isFALSE(all.equal(sum(c(0.3,p)),1))) stop("'p0' and 'p' must sum to 1")

  if(!is.numeric(vcv)) stop("'vcv' must be numeric")
  if(is.vector(vcv)){
    if(length(vcv)!=n) stop("number of variance elements in 'vcv' must match number of components")
    vcv <- array(rep(c(1,0,0,1),n)*rep(vcv^2,each=4),dim=c(2,2,n))
  } else if(is.array(vcv)){
    if(!all(dim(vcv)==c(2,2,n))) stop("'vcv' must be an array of 2x2 matrices with layers matching the number of components")
    for(i in 1:n){
      if((!isSymmetric(vcv[,,i]))||(det(vcv[,,i])<=0)) stop(paste("matrix",i,"in 'vcv' is invalid -- each must be symmetric and positive-definite"))
    }
  }

  w <- as.mask(window,dimyx=resolution)
  x <- w$xcol
  cellarea <- w$xstep*w$ystep
  y <- w$yrow
  xy <- expand.grid(x,y)
  xyinside <- inside.owin(x=xy[,1],y=xy[,2],w=w) # flags those coordinates of the pixel centroids that are actually inside W
  xyin <- xy[xyinside,]

  # Sets f to be the first density function, appropriately scaled
  f <- dmvnorm(xy,mean=mean[,1],sigma=vcv[,,1])
  scale1 <- sum(f*xyinside*cellarea)

  if(scale1<0.01 && !inside.owin(mean[1,1],mean[2,1],w)) warning("Component 1 may be out of range")
  f <- p[1] * f / scale1 # Scales f so this distribution integrates to p[1]

  # Adds the other density functions
  if(length(mean[1,])>1){
    for(i in 2:length(mean[1,])) {
      fAdd <- dmvnorm(xy,mean=mean[,i],sigma=vcv[,,i])
      scale <- sum(fAdd*xyinside*cellarea)
      if(scale<0.01) warning(paste("Component",i,"may be out of range"))
      fAdd <- p[i] * fAdd / scale
      f <- f + fAdd
    }
  }

  #Adds uniform density
  f <- f + p0/area.owin(window)
  f.im <- im(t(matrix(f,nrow=resolution)),xcol=x,yrow=y)
  f.im <- f.im[w,drop=FALSE]*int

  return(f.im)
}
