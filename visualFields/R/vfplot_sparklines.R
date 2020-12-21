vfplot_sparklines <- function( vf, ylim = c( -5, 35 ), collin = NULL, ... ) {

  texteval <- "vfsettings$locini"
  locini   <- eval( parse( text = texteval ) )
  # get x and y locations
  texteval   <- paste( vf$tperimetry[1], "locmap$",  vf$tpattern[1], sep = "" )
  patternMap <- eval( parse( text = texteval ) )
  patternMap <- patternMap[,c( "xod", "yod" )]
  if( vf$seye[1] == "OS" ) patternMap$xod <- -patternMap$xod

  # separate x vals from y vals and remove blind spot
  lt      <- as.POSIXlt( vf$tdate )
  mon     <- lt$year * 12 + lt$mon
  yeardif <- ( mon - mon[1] ) / 12
  vals    <- vf[,locini:ncol( vf )]

  # get blind spot
  evaltxt <- paste( "vfsettings$", vf$tpattern[1], "$bs", sep = "")
  bs <- eval( parse( text = evaltxt ) )
  if( !is.na( bs[1] ) ) {
    patternMap <- patternMap[-bs,]
    vals       <- vals[,-bs]
  }
  if( is.null( collin ) ) collin <- rep( "gray10", length( patternMap$xod ) )
  
  oparfig <- par()$fig # save graphics state before setting to individual "pointwise" regions
  
  # plot regions in NDC coordinates NB: need to modify for non-24-2 patterns 
  locs <- cbind( grconvertX( patternMap$xod - 2, to = "ndc" ), grconvertX( patternMap$xod + 2, to = "ndc" ),
                 grconvertY( patternMap$yod - 3, to = "ndc" ), grconvertY( patternMap$yod + 3, to = "ndc" ) )
  
  #  set a new plot region for each of the sparklines
  for( i in 1:nrow( locs ) )  {
    par( fig = locs[i,], mar = c( 0, 0, 0, 0 ), mgp = c( 0, 0, 0 ), new = TRUE )
    plot( yeardif, vals[,i], ylim = ylim, type = "l", axes = FALSE, col = collin[i] )
  }
  par ( fig = oparfig ) # reset the graphics states 
}