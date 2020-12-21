vftessellation <- function( locmap, dist = 3 ) {

  locmap <- locmap[,c( "xod", "yod" )]

  dxy    <- deldir( locmap$xod, locmap$yod )
  plocmap <- ppp( locmap$xod, locmap$yod,
                  c( min( locmap$xod ), max( locmap$xod ) ),
                  c( min( locmap$yod ), max( locmap$yod ) ) )
  del <- delaunay( plocmap )
  x   <- del$window$bdry[[1]]$x
  y   <- del$window$bdry[[1]]$y
  
  len <- sqrt( x^2 + y^2 ) + dist 
  ang <- atan2( y, x )

  CP <- list( x = len * cos( ang ), y = len * sin( ang ) )

  return( list( dxy, CP ) )
}