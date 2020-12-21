stintegral <- function(x){
  if(inherits(x,"stim")){
    a <- x$a
    xx <- x$xcol
    yy <- x$yrow
    tt <- x$tlay
  } else if(inherits(x,"stden")){
    a <- abind(lapply(x$z,as.matrix),along=3)
    xx <- x$z[[1]]$xcol
    yy <- x$z[[1]]$yrow
    tt <- x$tgrid
  } else {
    stop("'x' must be of class 'stim' or 'stden'")
  }
  return(sum(a*(xx[2]-xx[1])*(yy[2]-yy[1])*(tt[2]-tt[1]),na.rm=TRUE))
}
