psi2oct <- function( psi, diam = 12 ) {
  radius <- diam / 2
  r <- seq( 4, 10, 0.001 )
  oct <- NULL
  for( i in 1:length( psi ) ) {
    fb     <- fiberpathpsi( psi[i] )
    xy     <- jpolar2cart( data.frame( r = r, phi = fb( r ) ) ) - c( 15, 2 )
    xy     <- xy[which.min( abs( xy$x^2 + xy$y^2 - radius^2 ) ),]
    oct2   <- atan2( xy$y, -xy$x )
    oct2[oct2 < 0] <- 2 * pi + oct2[oct2 < 0]
    oct[i] <- oct2
  }
  return( 180 / pi * oct )
}