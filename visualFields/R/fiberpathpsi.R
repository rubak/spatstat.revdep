fiberpathpsi <- function( psi0, r0 = 4 ) {
  # Returns a function describing the expected fiber path given an angle of
  # incidence on the blind spot. Implemented after Jansonius_et_al_VR_2009, but see
  # erratum at Jansonius_etal_VR_2010_erratum. Includes the nasal sectors from
  # Jansonius_et_al_EER_2012
  if( psi0 >= 0 & psi0 < 60 ) {
    b <- 0.00083 * psi0^2 + 0.020 * psi0 - 2.65
    c <- 1.9 + 1.4 * tanh( ( psi0 - 121 ) / 14 )
  }
  if( psi0 >= 60 & psi0 <= 180 ) {
    b <- exp( -1.9 + 3.9 * tanh( - ( psi0 - 121 ) / 14 ) )
    c <- 1.9 + 1.4 * tanh( ( psi0 - 121 ) / 14 )
  }
  if( psi0 > -180 & psi0 <= -60 ) {
    b <- -exp( 0.7 + 1.5 * tanh( - ( - psi0 - 90 ) / 25 ) )
    c <- 1.0 + 0.5 * tanh( ( -psi0 - 90 ) / 25 )
  }
  if( psi0 > -60 & psi0 < 0 ) {
    b <- 0.00083 * psi0^2 + 0.020 * psi0 - 2.65
    c <- 1.0 + 0.5 * tanh( ( -psi0 - 90 ) / 25 )
  }
  fbfun <- function( r ) {
    return( psi0 + b * ( r - r0 )^c )
  }
  return( fbfun )
}