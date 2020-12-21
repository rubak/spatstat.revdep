#
# linalg.R
#
# Miscellaneous linear algebra
#
#  $Revision: 1.1 $ $Date: 2013/07/16 05:45:42 $ 


# square root of pos def matrix
sqrtmat <- function(M) {
  s <- svd(M)
  d <- s$d
  n <- length(d)
  dsd <- diag(sqrt(d), n, n)
  Y <- s$u %*% dsd %*% t(s$v)
  return(Y)
}
# inverse square root of pos def matrix
invsqrtmat <- function(M) {
  s <- svd(M)
  d <- s$d
  n <- length(d)
  isd <- diag(1/sqrt(d), n, n)
  Y <- s$u %*% isd %*% t(s$v)
  return(Y)
}
