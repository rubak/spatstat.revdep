vfplot_plr <- function( sl, pval, vfinfo, newWindow = FALSE,
                        xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL,
                        colorMapType = "pval", colorScale = NULL,
                        txtfont = "sans", pointsize = 10, width = 6,
                        showaxis = FALSE, colaxis = "black" ) {

  # construct patternmap
  evaltxt    <- paste( vfinfo$tperimetry, "locmap$", vfinfo$tpattern, sep = "" )
  patternMap <- eval( parse( text = evaltxt ) )
  patternMap <- patternMap[, c( "xod", "yod" )]
  # get normative values
  texteval <- "vfenv$nv"
  nv       <- eval( parse( text = texteval ) )
  # get blind spot
  evaltxt <- paste( "vfsettings$", vfinfo$tpattern[1], "$bs", sep = "")
  bs <- eval( parse( text = evaltxt ) )
  if( !is.na( bs[1] ) ) {
    bsloc      <- patternMap[bs,]
    patternMap <- patternMap[-bs,]
  } else {
    bsloc <- data.frame( xod = c( 15, 15 ), yod = c( 3, -3 ) )
  }

  # types of color map and ring map
  if( is.null( colorMapType ) ) stop( "colorMapType must be 'slope', 'pval', or 'blind'" )
  if( colorMapType != "pval" & colorMapType != "slope" & colorMapType != "blind" ) stop( "wrong colorMapType. Must be 'slope', 'pval', or 'blind'" )

  # get the conventional color scale
  if( colorMapType == "pval" ) {
    if( is.null( colorScale ) ) {
      colorScale  <- nv$pmapsettings
    }
    pval  <- 100 * pval
    pvalc <- rep( 100, length( pval ) )
    pvalc[which( pval <= colorScale$cutoffs[1] )] <- colorScale$cutoffs[1]
    for( i in 2:( length( colorScale$cutoffs ) - 1 ) ) pvalc[which( pval > colorScale$cutoffs[i - 1] & pval <= colorScale$cutoffs[i] )] <- colorScale$cutoffs[i]
    valForMapping <- pvalc
  }
  # inform the color scale for slopes
  if( colorMapType == "slope" ) {
    if( is.null( colorScale ) ) {
      colorScale         <- NULL
      colorScale$cutoffs <- c( -1.5, -1.0, -0.5, 0.5, 1 )
      colorScale$red     <- c( 0.8914642, 0.9999847, 0.9999847, 0.9742432, 0.0000000 )
      colorScale$green   <- c( 0.0000000, 0.5706177, 0.9041748, 0.9355011, 0.9999847 )   
      colorScale$blue    <- c( 0.1622925, 0.1513214, 0.0000000, 0.9213409, 0.9999847 )
      colorScale         <- as.data.frame( colorScale )
    }
    # map slope values to corresponding categories defined by colorScale$cutoffs
    slc <- NULL
    slc[c( 1:length( sl ) )] <- NA
    slc[which( sl <= colorScale$cutoffs[1] )] <- colorScale$cutoffs[1]
    slc[which( sl > colorScale$cutoffs[length( colorScale$cutoffs )] )] <- colorScale$cutoffs[length( colorScale$cutoffs )]
    for( i in 2:length( colorScale$cutoffs ) ) {
      slc[which( sl > colorScale$cutoffs[i-1] & sl <= colorScale$cutoffs[i] )] <- colorScale$cutoffs[i]
    }
    valForMapping <- slc
  }
  
  plotColor <- vfcolormap( valForMapping, mapval = colorScale )

  # if not imposed, calculate limits of plot expand by 5% each axis
  xrange <- max( patternMap$xod ) - min( patternMap$xod )
  yrange <- max( patternMap$yod ) - min( patternMap$yod )
  if( is.null( xmin ) ) xmin <- min( patternMap$xod ) - 0.025 * xrange
  if( is.null( xmax ) ) xmax <- max( patternMap$xod ) + 0.025 * xrange
  if( is.null( ymin ) ) ymin <- min( patternMap$yod ) - 0.025 * yrange
  if( is.null( ymax ) ) ymax <- max( patternMap$yod ) + 0.025 * yrange

  # add void points if not in the test of locations
  if( length( bs ) > 0 ) {
    patternMap <- rbind( patternMap, bsloc )
    plotColor  <- rbind( plotColor, data.frame( red = c( 0.5, 0.5 ), green = c( 0.5, 0.5 ), blue =  c( 0.5, 0.5 ) ) )
    sl         <- c( sl, rep( NA, 2 ) )
  }

  if( vfinfo$seye == "OS" ) {
    xmin2 <- xmin
    xmin  <- -xmax
    xmax  <- -xmin2
    patternMap$xod <- -patternMap$xod
  }

  # get the Voronoi tesselation tiles tiles
  vftess  <- vftessellation( patternMap, dist = 3 )
  vftiles <- tile.list( vftess[[1]] )
  vfhull  <- vftess[[2]]

  # create a new window and plot data in it
  # window rescale is set to fixed to ensure re-sizing window doesn't re-size the plot
  if( newWindow ) {
    height <- width * ( ymax - ymin ) / ( xmax - xmin )
    device <- options( "device" )
    if( .Platform$OS.type == "unix" ) {
      if( Sys.info()["sysname"] == "Darwin" ) {
        options( device = "quartz" )
        dev.new( width = width, height = height, dpi = 85 )
      } else {
        options( device = "x11" )
        dev.new( width = width, height = height )
      }
    } else{
      options( device = "windows" )
      dev.new( width = width, height = height, rescale = "fixed" )
    }
    options( device = device )
  }

# slope is in dB per 10 years
  sl <- round( 10 * sl, 1 )

  vfplotloc( sl, patternMap = patternMap,
             vftiles = vftiles, vfhull = vfhull, loccol = plotColor,
             xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
             txtfont = txtfont, pointsize = pointsize,
             showaxis = showaxis, colaxis = colaxis )
}