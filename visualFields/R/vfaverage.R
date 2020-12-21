vfaverage <- function( vf ) {

  # check input
  if( length( unique( vf$tperimetry ) ) > 1 |
      length( unique( vf$tpattern   ) ) > 1 |
      length( unique( vf$talgorithm ) ) > 1 |
      length( unique( vf$id ) ) > 1         |
      length( unique( vf$seye ) ) > 1 ) {
    stop( "all visual fields should belong to the same subject and eye tested with the same perimeter and algorithm on the same locations" )
  }

  # average values
  texteval <- "vfsettings$locini"
  locini   <- eval( parse( text = texteval ) )

  # get info data
  vfmean           <- vf[1,]
  vfmean$tdate     <- mean( vf$tdate )
  vfmean$ttime     <- NA
  vfmean$sage      <- mean( vf$sage )
  vfmean$sbsx      <- mean( vf$sbsx )
  vfmean$sbsy      <- mean( vf$sbsy )
  vfmean$sfp       <- mean( vf$sfp )
  vfmean$sfn       <- mean( vf$sfn )
  vfmean$sfl       <- mean( vf$sfl )
  vfmean$sduration <- NA
  vfmean$spause    <- NA

  # average per location
  vfmean[,locini:ncol( vf )] <- colMeans( vf[,locini:ncol( vf )] )

  return( vfmean )
  
}