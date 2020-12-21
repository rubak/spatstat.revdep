unify.owin <- function(W){
  if(!is.owin(W)) stop("'W' must be of spatstat class 'owin'")
  xr <- W$xrange
  yr <- W$yrange
  dx <- diff(xr)
  dy <- diff(yr)
  return(affine(W,matrix(c(1/dx,0,0,1/dy),2,2),c(-xr[1]/dx,-yr[1]/dy)))
}
