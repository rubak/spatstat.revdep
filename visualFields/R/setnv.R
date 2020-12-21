setnv <- function( nvtxt = "nvsapdefault" ) {

  nvinternal        <- eval( parse( text = nvtxt ) ) 
  nvinternal$nvname <- nvtxt
  texteval          <- "assign( \"nv\", nvinternal, envir = vfenv )"

  eval( parse( text = texteval ) )

}
