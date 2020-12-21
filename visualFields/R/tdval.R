tdval <- function( vf ) {

  texteval <- "vfsettings$locini"
  locini   <- eval( parse( text = texteval ) )

# get normal age-corrected values
  td <- vf
  for( i in 1:nrow( td ) ) {
# get how many locations we need to look at
    texteval <- paste( "vfsettings$", td$tpattern[i], "$locnum", sep = "" )
    locnum <- eval( parse( text = texteval ) )
# get blind-spot position
    texteval <- paste( "vfsettings$", td$tpattern[i], "$bs", sep = "" )
    bspos <- eval( parse( text = texteval ) )
# get the norm data and calculate normal age-corrected sensitivities
    texteval <- paste( "vfenv$nv$", td$tpattern[i], "_", td$talgorithm[i], "$agelm", sep = "" )
    agelm <- eval( parse( text = texteval ) )
# calculate total-deviation values
    td[i,locini:( locini - 1 + locnum )] <- td[i,locini:( locini - 1 + locnum )] - ( agelm$intercept + agelm$slope * td$sage[i] )
  }

# Fill with NA the blind-spot sensitivities
  if( all( !is.na( bspos[1] ) ) ) td[,bspos + locini - 1] <- NA 

  return( td )
}