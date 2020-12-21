vfplot <- function( vf, plotType = "vf",
                    xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL,
                    notSeenAsBlack = TRUE, newWindow = FALSE,
                    txtfont = "sans", pointsize = 10, width = 6,
                    showaxis = FALSE, colaxis = "black" ) {

  # check that vf has only 1 entry
  if( nrow( vf ) > 1 ) {
    stop("Error! vf cannot have more than 1 rows")
  }

  texteval <- "vfsettings$locini"
  locini   <- eval( parse( text = texteval ) )
  
  # construct the pattern string based on the pattern type
  evaltxt <- paste( "vfsettings$", vf$tpattern, "$locnum", sep = "" )
  loc_num <- eval( parse( text = evaltxt ) )

  # construct patternmap
  evaltxt    <- paste( vf$tperimetry, "locmap$", vf$tpattern, sep = "" )
  patternMap <- eval( parse( text=evaltxt ) )
  patternMap <- patternMap[,c( "xod", "yod" )]

  # get blind spot position
  evaltxt <- paste( "vfsettings$", vf$tpattern, "$bs", sep = "" )
  bs      <- eval( parse( text = evaltxt ) )

  # special locations in the visual field: BS locations
  bsxy <- data.frame( xod = c( 15, 15), yod = c( 3, -3 ) )

  # if not imposed, calculate limits of plot
  # expand by 5% each axis
  xrange <- max( patternMap$xod ) - min( patternMap$xod )
  yrange <- max( patternMap$yod ) - min( patternMap$yod )
  if( is.null( xmin ) ) xmin <- min( patternMap$xod ) - 0.025 * xrange
  if( is.null( xmax ) ) xmax <- max( patternMap$xod ) + 0.025 * xrange
  if( is.null( ymin ) ) ymin <- min( patternMap$yod ) - 0.025 * yrange
  if( is.null( ymax ) ) ymax <- max( patternMap$yod ) + 0.025 * yrange

# if plot type is "td", "pd", or "pdghr" then calculate corresponding deviation maps
  if( plotType == "td" ) {
    dev <- tdval( vf )
    vfp <- tdpmap( dev )
  }  
  if( plotType == "pd" ) {
    dev <- pdval( tdval( vf ) )
    vfp <- pdpmap( dev )
  }
  if( plotType == "pdghr" ) {
    dev <- pdvalghr( tdval( vf ) )
    vfp <- pdpmapghr( dev )
  }

# get colors for TD, PD maps or grayscale with VF plots
  if( plotType == "vf" ) {
    plotColor <- vfgrayscale( vf[,locini:( locini + loc_num - 1 )], age = vf$sage, pattern = vf$tpattern, algorithm = vf$talgorithm )
    vals <- as.character( round( vf[,locini:( locini + loc_num - 1 )] ) )
    vals[which( vf[,locini:( locini + loc_num - 1 )] < 0 )] = "<0"
  }  else {
    plotColor <- vfcolormap( as.numeric( vfp[,locini:( locini + loc_num - 1 )] ) )
    vals <- as.character( round( as.numeric( dev[,locini:( locini + loc_num - 1 )] ) ) )
    if( notSeenAsBlack ) {
      idxblack <- which( vf[locini:( locini + loc_num - 1 )] <= 0)
      if( length( idxblack ) > 0 ) plotColor[idxblack,] <- 0
    }
  }

  # add blind spot locations if required
  if( is.na( bs[1] ) & 
      ( ( xmin < bsxy$x[1] & xmax > bsxy$x[1] & ymin < bsxy$y[1] & ymax > bsxy$y[1] ) |
        ( xmin < bsxy$x[2] & xmax > bsxy$x[2] & ymin < bsxy$y[2] & ymax > bsxy$y[2] ) ) ) {
    vals       <- c( vals, c( "", "" ) )
    patternMap <- rbind( patternMap, bsxy )
    if( plotType == "vf" ) col <- c( 0, 0, 0 )
    else                   col <- c( 0.5, 0.5, 0.5 )
    plotColor  <- rbind( plotColor, col )
    plotColor  <- rbind( plotColor, col )
  }

  if( vf$seye == "OS" ) {
    xmin2 <- xmin
    xmin  <- -xmax
    xmax  <- -xmin2
    patternMap$xod <- -patternMap$xod
  }
  # get the Voronoi tesselation tiles tiles
  vftess  <- vftessellation( patternMap, dist = 3 )
  vftiles <- tile.list( vftess[[1]] )
  vfhull  <- vftess[[2]]

  # remove BS values and locations from patternMap
  if( !is.na( bs[1] ) & plotType != "vf" ) {
    vals[bs] <- ""
    plotColor[bs,] <- c( 0.5, 0.5, 0.5 )
  }
# create a new window and plot data in it
# window rescale is set to fixed to ensure re-sizing window doesn't re-size the plot
  height <- width * ( ymax - ymin ) / ( xmax - xmin )
  if( newWindow ) {
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
  vfplotloc( vals, patternMap = patternMap, loccol = plotColor,
             vftiles = vftiles, vfhull = vfhull,
             xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
             txtfont = txtfont, pointsize = pointsize,
             showaxis = showaxis, colaxis = colaxis )
}