cart2jpolar <- function( xy ) {
  x  <- xy[,1]
  y  <- xy[,2]
  xp <- x - 15
  yp <- y
  yp[x > 0] <- y[x > 0] - 2 * ( x[x > 0] / 15 )^2

  r   <- sqrt( xp^2 + yp^2 )
  psi <- atan2( yp, xp )
  return( data.frame( r = r, psi = 180 / pi * psi ) )
}