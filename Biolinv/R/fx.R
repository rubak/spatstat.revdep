#' One dimentional dispersal kernel function.
#'
#' Computes user-defined one dimensional dispersal kernels.
#'
#' @param x vector of distances.
#' @param a Alpha value.
#' @param c C value.
#'
#' @return numeric vector.
#'
#' @author Luca Butikofer
#'
#' @export
#'
#' @examples
#' plot(1:1000, fx(x= 1:1000, a= 300, c=2), type='l',
#'  ylab= 'Probablitiy', xlab='Distance',
#'  main= 'One Dimensional Dispersal Kernel \n Alpha= 300; C= 2')


fx<-function(x,a,c){
  y<-(c/(2*a*gamma(1/c)))*exp(-(x/a)^c)
  return(y)
}
