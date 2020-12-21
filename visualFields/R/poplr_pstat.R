poplr_pstat <- function( vf, porder, sltest = NULL ) {
  ##############
  # input checks
  ##############
  # check that all rows in vf belong to the same subject, the same test, the same perimetry
  # testing and the same algorithm
  if( length( unique( vf$tperimetry ) ) > 1 |
      length( unique( vf$tpattern   ) ) > 1 |
      length( unique( vf$talgorithm ) ) > 1 |
      length( unique( vf$id ) ) > 1         |
      length( unique( vf$seye ) ) > 1 ) {
    stop( "all visual fields should belong to the same subject tested with the same perimeter and algorithm on the same locations" )
  }
  if( nrow( porder ) > 1000000 ) stop( "please don't! Don't use more than a million permutations!" )
  texteval <- "vfsettings$locini"
  locini   <- eval( parse( text = texteval ) )

  # extract years and location values from vf and delete blind spots
  years <- as.numeric( vf$tdate - vf$tdate[1] ) / 365.25 # it should be difference in years from basline date
  evaltxt <- paste( "vfsettings$", vf$tpattern[1], "$bs", sep = "" )
  bs <- eval( parse( text = evaltxt ) )
  vf <- vf[,locini:ncol( vf )]
  if ( !is.na( bs[1] ) ) vf <- vf[,-bs]
  vf <- as.matrix( vf )
  # number of permutations, locations, and tests
  nperm <- nrow( porder )
  nloc  <- ncol( vf )
  ntest <- nrow( vf )
  # init
  precision <- 1e-6
  res       <- NULL
  res$sl    <- matrix( c( NA ), nrow = nperm, ncol = nloc )
  res$int   <- matrix( c( NA ), nrow = nperm, ncol = nloc )
  res$se    <- matrix( c( NA ), nrow = nperm, ncol = nloc )
  # add defaults for slope hypothesis tests when slr analysis is to be performed
  if( is.null( sltest ) ) sltest <- rep( c( 0 ), nloc )
  # get the locations for which sensitivity did not change
  invariantloc <- as.numeric( which( colSds( vf ) <= precision ) )
  # point-wise linear regression over time permutation-invarian values
  syears  <- sum( years )
  myears  <- mean( years )
  ssyears <- ( ntest - 1 ) * var( years )
  kvyears <- ( ntest - 2 ) * ssyears
  mvf   <- c( colMeans( vf ) )
  ssvf  <- c( ( ntest - 1 ) * colVars( vf ) )
  # calculate regression slopes, intercepts, and slope standard errors per location
  for( loc in 1:nloc ) {
    res$sl[,loc]  <- ( matrix( years[porder], nrow( porder ), ncol( porder ) ) %*% vf[,loc]
                       - syears * mean( vf[,loc] ) ) / ssyears
    res$int[,loc] <- rep( mvf[loc], nperm ) - myears * res$sl[,loc]
    varslope <- ( rep( ssvf[loc], nperm ) - ssyears * res$sl[,loc]^2 ) / kvyears
    varslope[which( varslope < 0 )] <- 0
    res$se[,loc] <- sqrt( varslope )
  }
  # locations with non-changing series in sensitivity: slope is zero,
  # intercept is not defined, and standard error is nominally very small
  res$sl[,invariantloc]  <- 0
  res$int[,invariantloc] <- vf[1,invariantloc]
  res$se[,invariantloc]  <- precision
  # test sensitivity slope lower than specified slope
  res$locpvals <- pt( ( res$sl - t( matrix( rep( sltest, nloc * nperm ), nloc, nperm ) ) ) / res$se, ntest - 2 )
  return( res )
}