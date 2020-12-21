riw <- function(a,A){
  return(solve(rw(a,solve(A))))
}

rw <- function(b,B){
  p <- nrow(B)
  Z <- matrix(0,p,p)
  diag(Z) <- sqrt(rchisq(p,b:(b-p+1)))
  pseq <- 1:(p-1)
  Z[rep(p*pseq,pseq)+unlist(lapply(pseq,seq))] <- rnorm(p*(p-1)/2)
  return(crossprod(Z%*%chol(B)))
}
