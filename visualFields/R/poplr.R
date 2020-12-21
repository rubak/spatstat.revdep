poplr <- function( vf, nperm = 5000, sltest = NULL, truncVal = 1 ) {
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
  if( nperm < 5 ) stop( "number of permutations lower than 5" )
  if( nperm > 1000000 ) stop( "please don't! Don't use more than a million permutations!" )
  # truncation must be between zero and one
  if( truncVal <= 0 | truncVal > 1 ) stop("truncation must be between 0 and 1")

  texteval <- "vfsettings$locini"
  locini   <- eval( parse( text = texteval ) )

  # permutation matrix
  porder <- make.permSpace( c( 1:nrow( vf ) ), nperm, return.permIDs = TRUE )$permID
  porder <- rbind( c( 1:nrow( vf ) ), porder )
  # get the p-value statitics of the permutation analysis ...
  pstat <- poplr_pstat( vf, porder = porder, sltest = sltest )
  # ... and the actual analysis
  cstat <- poplr_cstat( pstat$locpvals, truncVal = truncVal )
  res <- NULL
  # get last VF in res$vfdata
  res$vfdata      <- vf[nrow( vf ),]
  # get and remove blind spot
  evaltxt <- paste("vfsettings$", vf$tpattern[1], "$bs", sep = "")
  bs <- eval( parse( text = evaltxt ) ) + locini - 1
  if ( !is.na( bs[1] ) ) res$vfdata <- res$vfdata[-bs]
  res$nvisits  <- nrow( vf )
  res$sltest   <- sltest
  res$nperm    <- nperm
  res$sl       <- pstat$sl[1,]
  res$int      <- pstat$int[1,]
  res$se       <- pstat$se[1,]
  res$locpvals <- pstat$locpvals[1,]
  res$s        <- cstat$s
  res$sp       <- cstat$sp
  res$pval     <- cstat$pval
  res$sr       <- cstat$sr
  res$spr      <- cstat$spr
  res$pvalr    <- cstat$pvalr
  return( res )
}