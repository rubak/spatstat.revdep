gcloc2psi <- function( xy, r0 = 4 ) {
  # Returns the angle of incidence at the ONH of a vf locations
  xy[,2] <- -xy[,2]
  psi0   <- rep( NA, nrow( xy ) )
  rpsi   <- cart2jpolar( xy )
  for( i in 1:nrow( xy ) ) {
    fbpathinv <- function( psi ) {
      fb <- fiberpathpsi( psi, r0 )
      psiest <- fb( rpsi[i,1] )
      return( ( psiest - rpsi[i,2] )^2 )
    }
    psi0[i] <- optimize( fbpathinv, interval = c( -179.99, 180 ) )$minimum
  }
  return( psi0 )
}