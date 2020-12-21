rgmix <- function(N,window,v=4,S=NULL,extras=FALSE,...){
  if(!is.owin(window)) stop("'window' must be of spatstat class 'owin'")

  if(!is.numeric(N)) stop("'N' must be numeric")
  N <- round(N[1])
  if(N<=0) stop("'N' must be positive")

  v <- round(v[1])
  if(v<4) stop("'v' must be at least 4")

  if(is.null(S)){
    wx <- window$xrange
    wy <- window$yrange
    S <- matrix(c(((wx[2]-wx[1])/12)^2,0,0,((wy[2]-wy[1])/12)^2),nrow=2)
  }

  if(!is.matrix(S)) stop("'S' must be a 2x2 matrix")
  if(!all(dim(S)==2)) stop("'S' must be a 2x2 matrix")
  if(!isSymmetric(S)||det(S)<=0) stop("'S' must be symmetric and positive-definite")

  pts <- runifpoint(N,window)
  mn <- rbind(pts$x,pts$y)
  vcvarray <- array(0,dim=c(2,2,N))
  for(i in 1:N) vcvarray[,,i] <- riw(v,S)

  result <- sgmix(mn,vcvarray,window,...)
  if(extras) result <- list(f=result,mn=mn,vcv=vcvarray)
  return(result)
}
