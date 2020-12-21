rrpoint <- function(n,r,W=NULL,correction=1.1,maxpass=50){

  if(!inherits(r,"rrim")){
    if(!inherits(r,"rrs")){
      stop("'r' must be of class 'rrim' or 'rrs'")
    } else {
      rs <- r
      r <- list(f=r$f$z,g=r$g$z)
      W <- Window(rs$f$pp)
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

  fdat <- rimpoly(n[1],z=r$f,w=W,correction=correction,maxpass=maxpass)
  gdat <- rimpoly(n[2],z=r$g,w=W,correction=correction,maxpass=maxpass)

  return(list(cases=fdat,controls=gdat))
}
