vf2gcloc <- function( xy ) {
  # get function of VF location to GC soma location
  ecc <- eval( parse( text = "gcdisp$ecc - gcdisp$displ" ) )
  displ <- eval( parse( text = "gcdisp$displ" ) )
  # convert to polar coordinates
  r     <- sqrt( xy[,1]^2 + xy[,2]^2 )
  theta <- atan2( xy[,2], xy[,1] )
  # get soma GC eccentricity
  r <- r + approx( ecc, displ, r, rule = 2 )$y
  #convert back to cartesian coordinates
  xy[,1] <- r * cos( theta )
  xy[,2] <- r * sin( theta )

  return( xy )
}