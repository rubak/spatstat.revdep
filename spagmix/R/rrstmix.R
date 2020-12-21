rrstmix <- function(g, rhotspots, rsds, rweights, rbase = 1, log = TRUE, tlim = NULL, tres = NULL) {

  if(is.im(g)){
    if(is.null(tlim)||is.null(tres)) stop("temporal interval 'tlim' and resolution 'tres' must be supplied when 'g' is a spatial 'im'")
    if(!is.vector(tlim)) stop("'tlim' must be a vector")
    if(length(tlim)!=2) stop("'tlim' must be a vector of length 2")
    if(!is.numeric(tlim)) stop("'tlim' must be numeric")
    if(tlim[2]<=tlim[1]) stop("'tlim[2]' must be greater than 'tlim[1]'")
    if(tres<=1) stop("'tres' must be > 1")

    tgr <- seq(tlim[1],tlim[2],length=tres+1)
    tlay <- tgr[-1] - 0.5*diff(tlim)/tres
    ga <- array(NA,dim=c(length(g$xcol),length(g$yrow),tres))
    for(i in 1:tres) ga[,,i] <- t(as.matrix(g))
    vol <- (g$xcol[2]-g$xcol[1])*(g$yrow[2]-g$yrow[1])*(tlay[2]-tlay[1])
    ga <- ga/sum(ga*vol,na.rm=TRUE)
    gv <- solist()
    for(i in 1:tres) gv[[i]] <- im(t(ga[,,i]),xcol=g$xcol,yrow=g$yrow)
    g <- list(a=ga,v=gv,xcol=g$xcol,yrow=g$yrow,tlay=tlay,W=g$W)
    class(g) <- "stim"
  } else {
    if(!inherits(g,"stim")){
      if(!inherits(g,"stden")){
        stop("'g' must be of class 'stim', 'stden', or 'im'")
      } else {
        f <- list(v=g$z,tlay=g$tgrid,a=abind(lapply(g$z,as.matrix),along=3))
        g <- f
      }
    }
  }

  wn <- g$W

  if(!is.matrix(rhotspots)||!is.numeric(rhotspots)) stop("'rhotspots' must be a numeric matrix")
  if(nrow(rhotspots)!=3) stop("'rhotspots' must have 3 rows")
  nh <- ncol(rhotspots)

  if(!is.matrix(rsds)||!is.numeric(rsds)) stop("'rsds' must be a numeric matrix")
  if(nrow(rsds)!=3) stop("'rsds' must have 3 rows")
  if(any(rsds<=0)) stop("'rsds' must be strictly positive")
  ns <- ncol(rsds)
  # ns <- length(rsds)

  if(!is.vector(rweights)||!is.numeric(rweights)) stop("'rweights' must be a numeric vector")
  nw <- length(rweights)

  if(!all(nh==c(ns,nw))) stop("'rhotspots', 'rsds', and 'rweights' must all have the same number of components")
  n <- nh

  rbase <- rbase[1]
  if(!is.numeric(rbase)) stop("'rbase' must be numeric")

  xcol <- g$xcol
  yrow <- g$yrow
  tlay <- g$tlay
  xn <- length(xcol)
  yn <- length(yrow)
  tn <- length(tlay)
  vol <- (g$xcol[2]-g$xcol[1])*(g$yrow[2]-g$yrow[1])*(g$tlay[2]-g$tlay[1])
  g$a <- g$a/sum(g$a*vol,na.rm=TRUE) # ensures g integrates to 1
  g$v <- lapply(g$v,function(x) x/sum(g$a*vol,na.rm=TRUE))

  xyt <- expand.grid(xcol,yrow,tlay)
  rarr <- rbase
  for(i in 1:n){
    u <- (xyt[,1]-rhotspots[1,i])^2/rsds[1,i]^2 + (xyt[,2]-rhotspots[2,i])^2/rsds[2,i]^2 + (xyt[,3]-rhotspots[3,i])^2/rsds[3,i]^2
    rAdd <- rweights[i]*exp(-0.5*u)
    rarr <- rarr + rAdd
  }
  rarr <- array(rarr,dim=c(yn,xn,tn))
  rvol <- sum(rarr*g$a*vol,na.rm=TRUE)
  rarr <- rarr/rvol
  rarr[rarr<=0] <- .Machine$double.xmin
  farr <- rarr*g$a

  narem <- g$a
  narem[!is.na(narem)] <- 1

  rarr <- rarr*narem
  if(log) rarr <- log(rarr)

  rlist <- flist <- solist()
  for(i in 1:tn){
    rlist[[i]] <- im(t(rarr[,,i]),xcol=xcol,yrow=yrow)#[window, drop=FALSE]
    flist[[i]] <- im(t(farr[,,i]),xcol=xcol,yrow=yrow)#[window, drop=FALSE]
  }
  r.result <- list(a=rarr,v=rlist,xcol=xcol,yrow=yrow,tlay=tlay,W=wn)
  f.result <- list(a=farr,v=flist,xcol=xcol,yrow=yrow,tlay=tlay,W=wn)
  class(r.result) <- "stim"
  class(f.result) <- "stim"

  result <- list(f=f.result,g=g,r=r.result)
  class(result) <- "rrstim"

  return(result)
}
