#' Compute the boundary of disks
#' 
#' @param x point pattern
#' @param r radius of the disks
#' 
#' @export

arcs <- function(x, r){
  x <- internalise_pp(x)
  if(x$dim==3) stop("Area fraction function only for 2d.")
  resA <- .External("SGCS_getArcs_c",
                    x,
                    r,
                    PACKAGE="SGCS")
  A<-matrix(resA, ncol=3, byrow=T)
  A[,1] <- A[,1] + 1
  colnames(A) <- c("point", "angle1", "angle2")
  A
}