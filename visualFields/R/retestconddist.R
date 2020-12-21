retestconddist <- function( vf, nbase = 1, nfollow = 1, alpha = 0.1, typequantile = 7 ) {
  
  if( any( vf$tperimetry != vf$tperimetry[1] ) ) {
    stop( "The perimetry test in all tests in vf should be the same" )
  }
  if( any( vf$talgorithm != vf$talgorithm[1] ) ) {
    stop( "The algorithm in all tests in vf should be the same" )
  }
  if( any( vf$tpattern != vf$tpattern[1] ) ) {
    stop( "The pattern of presented locations in all tests in vf should be the same" )
  }
  if( nbase < 1 | nbase > 3 ) {
    stop( "number of baseline series between 1 and 3 please" )
  }
  if( nfollow < 1 | nfollow > 3 ) {
    stop( "number of followup series between 1 and 3 please" )
  }

  # first sort vfs
  vf <- vfsort( vf )
  # check if data is right
  idu <- as.data.frame( unique( cbind( vf$id, vf$seye ) ) )
  names( idu ) <- c( "id", "seye" )
  # get number of subject/eyes
  neyes <- nrow( idu )
  for( i in 1:neyes ) {
    idx <- which( vf$id == idu$id[i] & vf$seye == idu$seye[i] )
    idu$visits <- length( idx )
  }
  if( any( idu$visits != ( nbase + nfollow ) ) ) stop( "wrong number of visits for at least 1 subject for the analysis" )

  texteval <- "vfsettings$locini"
  locini   <- eval( parse( text = texteval ) )
  # get blind spot
  txteval <- paste( "vfsettings$", vf$tpattern[1], "$bs", sep = "" )
  bs <- eval( parse( text = txteval ) )
  # remove blind spot as needed
  if( !is.na( bs[1] ) ) {
    vf <- vf[,-c( locini + bs - 1 )]
  }
  # check all possible values
  retestcdist  <- NULL
  retestcdist$x <- sort( unique( c( as.matrix( vf[,locini:ncol( vf )] ) ) ) )
  retestcdist$y <- as.list( rep( NA, length( retestcdist$x ) ) )
  # construct all possible permutations from the series
  idxv <- seq( from = 1, to = ( nbase + nfollow ) )
  combs <- combinations( ( nbase + nfollow ), nbase )
  # start analysis for each possible combination
  for( i in 1:nrow( combs ) ) {
    # do the analysis for each eye
    for( j in 1:neyes ) {
      # first get the baseline data vs the followup data
      idx      <- which( vf$id == idu$id[j] & vf$seye == idu$seye[j] )
      dbbase   <- vf[idx[combs[i,]],locini:ncol( vf )]
      dbfollow <- as.numeric( colMeans( vf[idx[-combs[i,]],locini:ncol( vf )] ) )
      # analyse dbbase
      idxbase <- NULL
      for( k in 1:ncol( dbbase ) ) {
        if( all( dbbase[,k] == dbbase[1,k] ) ) idxbase <- c( idxbase, k )
      }
      # mount conditional distribution
      for( k in 1:length( idxbase ) ) {
        idxcp <- which( retestcdist$x == dbbase[1,idxbase[k]] )
        if( length( retestcdist$y[[idxcp]] ) == 1 && is.na( retestcdist$y[[idxcp]] ) ) {
          retestcdist$y[[idxcp]] <- dbfollow[idxbase[k]]
        } else {
          retestcdist$y[[idxcp]] <- c( retestcdist$y[[idxcp]], dbfollow[idxbase[k]] )
        }
      }
    }
  }
  # remove whatever has no data in it
  idx <- which( is.na( retestcdist$y ) )
  if( length( idx ) > 0 ) {
    retestcdist$x      <- retestcdist$x[-idx]
    retestcdist$y[idx] <- NULL
  }
  # sort the results and obtain percentiles at level alpha
  for( i in 1:length( retestcdist$y ) ) {
    retestcdist$y[[i]]  <- sort( retestcdist$y[[i]] )
    qq <- as.numeric( quantile( retestcdist$y[[i]], probs = c( alpha / 2, 1 - alpha / 2 ), type = typequantile ) )
    retestcdist$n[i]    <- length( retestcdist$y[[i]] )
    retestcdist$ylow[i] <- qq[1]
    retestcdist$yup[i]  <- qq[2]
  }
  return( retestcdist )
}