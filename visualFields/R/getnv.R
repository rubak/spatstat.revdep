getnv <- function() {

  texteval <- "vfenv$nv"
  return( eval( parse( text = texteval ) ) )

}
