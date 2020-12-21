jpolar2cart <- function( rpsi ) {
  r   <- rpsi[,1]
  psi <- rpsi[,2]
  psi <- ( pi / 180 * psi ) %% ( 2 * pi )
  x <- sqrt( r^2 / ( 1 + ( tan( psi ) )^2 ) )
  y <- abs( x * tan( psi ) )
  x[psi > ( pi / 2 ) & psi < ( 3 * pi / 2 )] <- -x[psi > ( pi / 2 ) & psi < ( 3 * pi / 2 )]
  y[psi > pi] <- -y[psi > pi]

  x        <- x + 15
  y[x > 0] <- y[x > 0] + 2 * ( x[x > 0] / 15 )^2
  return( data.frame( x = round( x, 6 ), y = round( y, 6 ) ) )
}