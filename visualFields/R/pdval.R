pdval <- function( td ) {

  texteval <- "vfsettings$locini"
  locini   <- eval( parse( text = texteval ) )

  pd <- td

  for( i in 1:nrow( pd ) ) {
# get how many locations we need to look at
    texteval <- paste( "vfsettings$", pd$tpattern[i], "$locnum", sep = "" )
    locnum <- eval( parse( text = texteval ) )
# get PD values from obtained gh
    pd[i,locini:( locini - 1 + locnum )] <- pd[i,locini:( locini - 1 + locnum )] - ghpostd( td[i,] )
  }

  return( pd )
}
