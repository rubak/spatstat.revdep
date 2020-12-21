rrstpoint <- function(n,r,W=NULL,correction=1.5,maxpass=50){

  if(!inherits(r,"rrstim")){
    if(!inherits(r,"rrst")){
      stop("'r' must be of class 'rrstim' or 'rrst'")
    } else {
      W <- Window(r$g$pp)
      if(!inherits(r$g,"stden")){
        xcol <- r$g$z$xcol
        yrow <- r$g$z$yrow
        tlay <- r$f$tgrid
        tres <- length(tlay)

        # print(tres)
        ga <- array(NA,dim=c(length(xcol),length(yrow),tres))
        # print(dim(ga))

        for(i in 1:tres) ga[,,i] <- t(as.matrix(g$z))
        vol <- (xcol[2]-xcol[1])*(yrow[2]-yrow[1])*(tlay[2]-tlay[1])
        ga <- ga/sum(ga*vol,na.rm=TRUE)
        gv <- solist()
        for(i in 1:tres) gv[[i]] <- im(t(ga[,,i]),xcol=xcol,yrow=yrow)
        g <- list(a=ga,v=gv,xcol=xcol,yrow=yrow,tlay=tlay,W=W)
        class(g) <- "stim"
        r$g <- g
      }
    }
  }

  if(!is.numeric(n)) stop("'n' must be numeric")
  nl <- length(n)
  if(nl!=1){
    if(nl==0) stop("'n' must contain at least one entry")
    n <- round(n[1:2])
  } else {
    n <- round(c(n,n))
  }
  if(any(n<=0)||any(is.infinite(n))) stop("entries of 'n' must be positive and finite")

  fdat <- rstpoint(n[1],f=r$f,W=W,correction=correction,maxpass=maxpass)
  gdat <- rstpoint(n[2],f=r$g,W=W,correction=correction,maxpass=maxpass)

  return(list(cases=fdat,controls=gdat))
}
