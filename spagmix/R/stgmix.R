
stgmix <- function(mean, vcv, window, tlim, p0=0, p=NULL, sres=128, tres=sres, int=1) {
  if(!is.owin(window)) stop("'window' must be of spatstat class 'owin'")
  w <- window

  if(!is.matrix(mean)) stop("'mean' must be a matrix")
  if(nrow(mean)!=3) stop("'mean' is of incorrect dimension")
  n <- ncol(mean)

  if(int<=0) stop("'int' must be positive")
  if(sres<=1||tres<=0) stop("'sres' and 'tres' must be >= 1")

  if(length(p0)>1) p0 <- p0[1]
  if(!is.numeric(p0)) stop("'p0' must be numeric")
  if(p0<0||p0>1) stop("'p0' must be in [0,1]")

  if(is.null(p)) p <- rep((1-p0)/n,n)
  if(!is.numeric(p)) stop("'p' must be numeric")
  if(length(p)!=n)
  if(any(p<0)||any(p>1)) stop("all elements of 'p' must be in [0,1]")
  if(sum(c(p,p0))!=1) stop("'p0' and 'p' must sum to 1")

  if(!is.numeric(vcv)) stop("'vcv' must be numeric")
  if(is.array(vcv)){
    if(!all(dim(vcv)==c(3,3,n))) stop("'vcv' must be an array of 3x3 matrices with layers matching the number of components")
    for(i in 1:n){
      if((!isSymmetric(vcv[,,i]))||(det(vcv[,,i])<=0)) stop(paste("matrix",i,"in 'vcv' is invalid -- each must be symmetric and positive-definite"))
    }
  } else stop("'vcv' must be an array")


  if(!is.vector(tlim)) stop("'tlim' must be a vector")
  if(length(tlim)!=2) stop("'tlim' must be a vector of length 2")
  if(!is.numeric(tlim)) stop("'tlim' must be numeric")
  if(tlim[2]<=tlim[1]) stop("'tlim[2]' must be greater than 'tlim[1]'")

  w <- as.mask(window,dimyx=sres)
  x <- w$xcol
  y <- w$yrow
  xy <- expand.grid(x,y)
  xyinside <- inside.owin(x=xy[,1],y=xy[,2],w=w) # flags those coordinates of the pixel centroids that are actually inside W
  xyin <- xy[xyinside,]
  tseq <- seq(tlim[1], tlim[2], length=tres+1)
  tstep <- (tlim[2] - tlim[1])/tres
  tt <- tseq[-(tres+1)] + 0.5*tstep #Voxel centroids along t-axis
  varea <- w$xstep*w$ystep*tstep
  xyt <- expand.grid(x,y,tt)


  f <- dmvnorm(xyt,mean=mean[,1],sigma=vcv[,,1])


  scale1 <- sum(f*xyinside*varea)

  if(scale1<0.01 && !inside.owin(mean[1,1],mean[2,1],w) && (mean[3,1]<tlim[1]||mean[3,1]>tlim[2])) warning("Component 1 may be out of range")
  f <- p[1] * f / scale1


  # Adds the other density functions
  if(n>1){
    for(i in 2:n){
      fAdd <- dmvnorm(xyt,mean=mean[,i],sigma=vcv[,,i])
      scale <- sum(fAdd*xyinside*varea)
      if (scale<0.01 && !inside.owin(mean[1,i],mean[2,i],w) && (mean[3,i]<tlim[1]||mean[3,i]>tlim[2])) warning(paste("Component",i,"may be out of range"))
      fAdd <- p[i] * fAdd / scale
      f <- f + fAdd
    }
  }

  #Adds uniform density
  volume <- sum(varea*sum(xyinside)*tres)
  f <- (f+p0/volume)*int

  narep <- rep(1,length(xyinside))
  narep[!xyinside] <- NA

  f.arr <- array(f*narep,dim=c(sres,sres,tres))
  f.list <- solist()
  for(i in 1:tres){
    f.list[[i]] <- im(t(f.arr[,,i]),xcol=x,yrow=y)
    f.list[[i]] <- f.list[[i]][w,drop=FALSE]
  }

  result <- list(a=f.arr,v=f.list,xcol=x,yrow=y,tlay=tt,W=window)
  class(result) <- "stim"
  return(result)
}
