rpoispoly <- function(z,w=NULL,correction=1.1,maxpass=50){
  if(is.im(z)){
    if(any(z<0)) stop("'z' must be a non-negatively-valued spatial intensity function")
  } else {
    stop("'z' must be an object of class 'im' (spatstat)")
  }
  return(rimpoly(rpois(1,integral(z)),z=z,w=w,correction=correction,maxpass=maxpass))
}
