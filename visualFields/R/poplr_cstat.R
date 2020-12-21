poplr_cstat <- function( pval, truncVal = 1 ) {
  ##############
  # input checks
  ##############
  # truncation must be between zero and one
  if( truncVal <= 0 | truncVal > 1 ) stop("truncation must be between 0 and 1")
  # init
  nperm     <- nrow( pval )
  nloc      <- ncol( pval )
  res       <- NULL
  res$s     <- NA
  res$sp    <- rep( NA, nperm )
  res$pval  <- NA
  res$sr    <- NA
  res$spr   <- rep( NA, nperm )
  res$pvalr <- NA

  # truncate p-values
  k <- matrix( rep( 1, nrow( pval ) * ncol( pval ) ), nrow( pval ), ncol( pval ) )
  k[pval > truncVal] <- 0
  kr <- matrix( rep( 1, nrow( pval ) * ncol( pval ) ), nrow( pval ), ncol( pval ) )
  kr[( 1 - pval ) > truncVal] <- 0

  # combine p-value test statistics
  # Fisher-class combination (product) of p-values with optional weigths
  res$sp  <- -rowSums( k  * log( pval ) )
  res$spr <- -rowSums( kr * log( 1 - pval ) )
  # observed and permutation test statistics
  res$s     <- res$sp[1]
  res$sr    <- res$spr[1]
  res$pval  <- 1 - rank( res$sp )[1]  / nperm
  res$pvalr <- 1 - rank( res$spr )[1] / nperm
  res$sp    <- res$sp[-1]
  res$spr   <- res$spr[-1]

  return( res )
}