#
#  First.R
#
# $Revision: 1.1 $  $Date: 2012/12/12 02:22:31 $
#

.onLoad <- function(...) { }

.onAttach <- function(libname, pkgname) {
  v <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                fields="Version")
  msg <- paste("spatstat.local", v)
  packageStartupMessage(msg)
  invisible(NULL)
}


