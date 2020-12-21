vfcolormap <- function( map, mapval = NULL ) {

  if( is.null( mapval ) ) {
    texteval <- "vfenv$nv$pmapsettings"
    mapval <- eval( parse( text = texteval ) )
  }

  rgbval                      <- NULL
  rgbval$red[1:length(map)]   <- c( NA )
  rgbval$green[1:length(map)] <- c( NA )
  rgbval$blue[1:length(map)]  <- c( NA )
  rgbval                      <- as.data.frame( rgbval )

  for( i in 1:length( mapval$cutoffs ) ) {
    idx <- which( map == mapval$cutoffs[i] )
    if( length( idx ) > 0 ) rgbval[idx,] <- mapval[i,2:ncol( mapval )]
  }
  return ( rgbval )
}
