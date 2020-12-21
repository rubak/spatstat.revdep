stkey <- function(start, stop, tlim, kf = NULL, tres = 64,
                  kftimes = NULL, fscales = NULL, window = NULL) {

  if(!is.im(start)||!is.im(stop)) stop("'start' and 'stop' must both be objects of spatstat class 'im'")

  if(!compatible.im(start,stop)) stop("pixel images not compatible")

  if(is.null(window)) window <- as.polygonal(Window(start))
  if(is.null(intersect.owin(Window(start),window,fatal=FALSE))) stop("'window' must overlap domain of supplied pixel images")

  start <- start / integral(start)
  stop <- stop / integral(stop)

  if(!is.numeric(tlim)||length(tlim)!=2) stop("'tlim' must be a numeric vector of length 2")

  if(tlim[2]<=tlim[1]) stop("'tlim' must be increasing")

  if(tres<2) stop("invalid 'tres'")

  nkey <- 0
  if(!is.null(kf)){
    if(!is.solist(kf)){
      stop("'kf' must be a spatstat 'solist' of pixel images")
    }
    nkey <- length(kf)
    comp <- sapply(kf,compatible.im,B=start)
    if(sum(comp)!=nkey) stop("pixel images in 'kf' must be compatible with 'start' and 'stop'")
  }
  kf <- lapply(kf,function(x) x/integral(x))

  ntot <- nkey+2
  if(nkey>0){
    if(is.null(kftimes)) kftimes <- seq(tlim[1],tlim[2],length=ntot)[-c(1,ntot)]
    if(length(kftimes)!=nkey) stop("length of 'kftimes' must match length of 'kf'")
    ftu <- unique(c(kftimes,tlim))
    if(length(ftu)!=(ntot)) stop("each spatial frame must correspond to a unique time given 'tlim' and 'ftimes'")
    if(any(kftimes>tlim[2])||any(kftimes<tlim[1])) stop("all 'kftimes' must fall within the interval given by 'tlim'")
    # get all times and spatial frames in increasing time order:
    ft <- sort(c(tlim,kftimes))
    kfto <- order(kftimes)
    f <- list()
    f[[1]] <- start
    for(i in 2:(ntot-1)) f[[i]] <- kf[kfto][[i-1]]
    f[[ntot]] <- stop
  } else {
    ft <- tlim
    f <- list(start,stop)
  }

  if(is.null(fscales)) fs <- fscales <- rep(1,ntot)
  if(any(fscales<0)) stop("'fscales' must be nonnegative")
  nfs <- length(fscales)
  if(nfs!=tres){
    if(nfs!=(ntot)) stop("an 'fscales' value must be specified for all spatial frames corresponding to c(start,kf,stop), *or* be of length 'tres'")
    if(nkey>0) fs <- c(fscales[1],fscales[kfto],fscales[ntot])
  } else {
    fs <- fscales
  }

  tseq <- seq(tlim[1],tlim[2],length=tres+1) #Edges of voxels on time axis
  tstep <- (tlim[2]-tlim[1])/tres
  tlay <- tseq[-(tres+1)] + 0.5*tstep #Voxel centroids on time axis

  arrayin <- abind(lapply(f,as.matrix),along=3)
  if(nfs==ntot){
    for(i in 1:nfs){
      arrayin[,,i] <- arrayin[,,i]*fs[i]
    }
    fmult <- rep(1,tres)
  } else {
    fmult <- fs
  }

  arrayout <- array(NA,dim=c(dim(start),tres))
  for(i in 1:dim(arrayout)[1]){
    for(j in 1:dim(arrayout)[2]){
      if(!is.na(start[i,j])){
        v <- arrayin[i,j,]
        arrayout[i,j,] <- approx(ft,v,xout=tlay)$y * fmult
      }
    }
  }

  volume <- sum(arrayout*start$xstep*start$ystep*tstep,na.rm=TRUE)
  arrayout <- arrayout/volume
  listout <- solist()
  for(i in 1:tres) {
    imout <- im(arrayout[,,i],xcol=start$xcol,yrow=start$yrow)
    listout[[i]] <- imout
  }

  result <- list(a=arrayout,v=listout,xcol=start$xcol,yrow=start$yrow,tlay=tlay,W=window)
  class(result) <- "stim"
  return(result)
}
