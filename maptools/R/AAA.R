.MAPTOOLS_CACHE <- new.env(FALSE, parent=globalenv())

#.onLoad <- function(lib, pkg) {
#    assign("gpclib", FALSE, envir=.MAPTOOLS_CACHE)
#}
register_s3_method <- function(pkg, generic, class, fun = NULL) {
  stopifnot(is.character(pkg), length(pkg) == 1L)
  stopifnot(is.character(generic), length(generic) == 1L)
  stopifnot(is.character(class), length(class) == 1L)

  if (is.null(fun)) {
    fun <- get(paste0(generic, ".", class), envir = parent.frame())
  } else {
    stopifnot(is.function(fun))
  }

  if (isNamespaceLoaded(pkg)) {
    registerS3method(generic, class, fun, envir = asNamespace(pkg))
  }

  # Always register hook in case package is later unloaded & reloaded
  setHook(
    packageEvent(pkg, "onLoad"),
    function(...) {
      registerS3method(generic, class, fun, envir = asNamespace(pkg))
    }
  )
}


.onLoad <- function(lib, pkg) {
    assign("gpclib", FALSE, envir=.MAPTOOLS_CACHE)
    rgeosI <- setRgeosStatus()
    if (getRversion() < "3.6.0") {
      register_s3_method("spatstat.geom", "as.im", "RasterLayer")
      register_s3_method("spatstat.geom", "as.im", "SpatialGridDataFrame")
      register_s3_method("spatstat.linnet", "as.linnet", "SpatialLines")
      register_s3_method("spatstat.geom", "as.owin", "SpatialGridDataFrame")
      register_s3_method("spatstat.geom", "as.owin", "SpatialPixelsDataFrame")
      register_s3_method("spatstat.geom", "as.owin", "SpatialPolygons")
      register_s3_method("spatstat.geom", "as.ppp", "SpatialPoints")
      register_s3_method("spatstat.geom", "as.ppp", "SpatialPointsDataFrame")
      register_s3_method("spatstat.geom", "as.psp", "Line")
      register_s3_method("spatstat.geom", "as.psp", "Lines")
      register_s3_method("spatstat.geom", "as.psp", "SpatialLines")
      register_s3_method("spatstat.geom", "as.psp", "SpatialLinesDataFrame")
    }
    invisible(NULL)
}

.onAttach <- function(lib, pkg) {
#    assign("gpclib", FALSE, envir=.MAPTOOLS_CACHE)
    Smess <- paste("Checking rgeos availability: ")
#    rgeosI <- setRgeosStatus()
    rgeosI <- rgeosStatus()
    Smess <- paste(Smess, rgeosI, "\n", sep="")
    if (!rgeosI) Smess <- paste(Smess, 
              "\tNote: when rgeos is not available, polygon geometry",
              "\tcomputations in maptools depend on gpclib,\n",
              "\twhich has a restricted licence. It is disabled by default;\n",
              "\tto enable gpclib, type gpclibPermit()\n")
    packageStartupMessage(Smess, appendLF = FALSE)
}

.onUnload <- function(libpath) {
    rm(.MAPTOOLS_CACHE)
}


