vflayout_progress <- function( vf, plotType, grp = 2, nperm = 5000,
                               colorMapType = "pval", colorScale = NULL,
                               filename = NULL,
                               pwidth = 8.27, pheight = 11.69, margin = 0.25,
                               showaxis = FALSE, colaxis = "black" ) {

  ##############
  # input checks
  ##############
  # check that all rows in vf belong to the same subject, the same test, the same perimetry
  # testing and the same algorithm, and the same eye
  if( length( unique( vf$tperimetry ) ) > 1 |
      length( unique( vf$tpattern   ) ) > 1 |
      length( unique( vf$talgorithm ) ) > 1 |
      length( unique( vf$id ) ) > 1         |
      length( unique( vf$seye ) ) > 1 ) {
    stop( "all visual fields should belong to the same subject and eye tested with the same perimeter and algorithm on the same locations" )
  }
  if( nrow( vf ) < 3 * grp ) stop( "the number of visual fields needs to be at least three times the number of visual fields to group for the analysis" )
  if( nrow( vf ) < 5 )       warning( "permutation analysis may not be very precise for less than 8 visual fields" )
  if( is.null( colorMapType) ) stop( "colorMapType must be 'slope', 'pval', or 'blind'" )
  if( colorMapType != "pval" & colorMapType != "slope" & colorMapType  != "blind" ) stop( "wrong colorMapType. Must be 'slope', 'pval', or 'blind'" )

  txtfont   <- "sans"
  pointsize <- 10

  # special locations in the visual field: BS locations
  bsxy   <- NULL
  bsxy$x <- c( 15, 15)
  bsxy$y <- c( 3, -3 )
  bsxy   <- as.data.frame( bsxy )
  
  # get settings for the pattern of test locations
  texteval <- "vfsettings$locini"
  locini   <- eval( parse( text = texteval ) )
  texteval <- "vfenv$nv"
  nv       <- eval( parse( text = texteval ) )
  texteval <- paste( "vfsettings$", vf$tpattern[1], sep = "" )
  settings <- eval( parse( text = texteval ) )
  # get x and y locations
  texteval   <- paste( vf$tperimetry[1], "locmap$",  vf$tpattern[1], sep = "" )
  patternMap <- eval( parse( text = texteval ) )
  patternMap <- patternMap[,c( "xod", "yod" )]

  # get the conventional color scale
  if( colorMapType == "pval" & is.null( colorScale ) ) {
    colorScale  <- nv$pmapsettings
  }
  if( colorMapType == "slope" & is.null( colorScale ) ) {
    colorScale         <- NULL
    colorScale$cutoffs <- c( -1.5, -1.0, -0.5, 0.5, 1 )
    colorScale$red     <- c( 0.8914642, 0.9999847, 0.9999847, 0.9742432, 0.0000000 )
    colorScale$green   <- c( 0.0000000, 0.5706177, 0.9041748, 0.9355011, 0.9999847 )   
    colorScale$blue    <- c( 0.1622925, 0.1513214, 0.0000000, 0.9213409, 0.9999847 )
    colorScale         <- as.data.frame( colorScale )
  }
  if( colorMapType == "blind" & is.null( colorScale ) ) {
    colorScale         <- NULL
    colorScale$cutoffs <- c( 5, 10, 15, 20 )
    colorScale$red     <- c( 0.8914642, 0.9999847, 0.9999847, 0.9742432 )
    colorScale$green   <- c( 0.0000000, 0.5706177, 0.9041748, 0.9355011 )   
    colorScale$blue    <- c( 0.1622925, 0.1513214, 0.0000000, 0.9213409 )
    colorScale         <- as.data.frame( colorScale )
  }

  ######################
  # analysis
  ######################
  # sort the data
  vf        <- vfsort( vf )
  # get global indices
  vfindices <- vfstats( vf )
  
  # get poplr analysis
  if( plotType == "vf" )    vals <- vf
  if( plotType == "td" )    vals <- tdval( vf )
  if( plotType == "pd" )    vals <- pdval( tdval( vf ) )
  if( plotType == "pdghr" ) vals <- pdvalghr( tdval( vf ) )
  
  pres <- poplr( vals, nperm = nperm )
  
  idxfirst    <- 1:grp
  idxmiddle   <- order( abs( vals$tdate - mean( c( vals$tdate[1], vals$tdate[nrow( vals )] ) ) ) )[1:grp]
  idxmiddle   <- idxmiddle[order( idxmiddle )]
  idxlast     <- ( nrow( vals ) - grp + 1 ):nrow( vals )
  tdatefirst  <- vals$tdate[idxfirst]
  tdatemiddle <- vals$tdate[idxmiddle]
  tdatelast   <- vals$tdate[idxlast]
  vffirst     <- vfaverage( vals[idxfirst,]  )
  vfmiddle    <- vfaverage( vals[idxmiddle,] )
  vflast      <- vfaverage( vals[idxlast,]   )

  if( plotType == "td" ) {
    vffirstp  <- tdpmap( vffirst  )
    vfmiddlep <- tdpmap( vfmiddle )
    vflastp   <- tdpmap( vflast   )
  }
  if( plotType == "pd" ) {
    vffirstp  <- pdpmap( vffirst  )
    vfmiddlep <- pdpmap( vfmiddle )
    vflastp   <- pdpmap( vflast   )
  }
  if( plotType == "pdghr" ) {
    vffirstp  <- pdpmapghr( vffirst  )
    vfmiddlep <- pdpmapghr( vfmiddle )
    vflastp   <- pdpmapghr( vflast   )
  }

  # open window wiht A4 page
  if( is.null( filename ) ) {
    device <- options( "device" )
    if( .Platform$OS.type == "unix" ) {
      if( Sys.info()["sysname"] == "Darwin" ) {
        options( device = "quartz" )
        dev.new( width = pwidth, height = pheight, dpi = 85 )
      } else {
        options( device = "x11" )
        dev.new( width = pwidth, height = pheight )
      }
    } else{
      options( device = "windows" )
      dev.new( width = pwidth, height = pheight, rescale = "fixed" )
    }
    options( device = device )
  } else {
    pdf( width = pwidth, height = pheight, file = filename )
  }

  # define the margins
  mwidth  <- pwidth  - 2 * margin
  mheight <- pheight - 2 * margin
  
  # create the layout of the printout
  printout <- createviewport( "printout", left = margin, top = margin, height = mheight, width = mwidth )

  xrange <- max( patternMap$xod ) - min( patternMap$xod )
  yrange <- max( patternMap$yod ) - min( patternMap$yod )
  xmin   <- min( patternMap$xod ) - 0.025 * xrange
  xmax   <- max( patternMap$xod ) + 0.025 * xrange
  ymin   <- min( patternMap$yod ) - 0.025 * yrange
  ymax   <- max( patternMap$yod ) + 0.025 * yrange

  if( plotType == "vf" ) {
    colorf   <- vfgrayscale( vffirst[,locini:ncol( vffirst )], age = vffirst$sage, pattern = vffirst$tpattern, algorithm = vffirst$talgorithm )
    vffirst  <- as.character( round( vffirst[,locini:ncol( vffirst )] ) )
    vffirst[which( vffirst < 0 )] = "<0"
    colorm   <- vfgrayscale( vfmiddle[,locini:ncol( vfmiddle )], age = vfmiddle$sage, pattern = vfmiddle$tpattern, algorithm = vfmiddle$talgorithm )
    vfmiddle <- as.character( round( vfmiddle[,locini:ncol( vfmiddle )] ) )
    vfmiddle[which( vfmiddle < 0 )] = "<0"
    colorl   <- vfgrayscale( vflast[,locini:ncol( vflast )], age = vflast$sage, pattern = vflast$tpattern, algorithm = vflast$talgorithm )
    vflast   <- as.character( round( vflast[,locini:ncol( vflast )] ) )
    vflast[which( vflast < 0 )] = "<0"
  } else {
    colorf   <- vfcolormap( as.numeric( vffirstp[,locini:ncol( vffirst )] ) )
    vffirst  <- as.character( round( as.numeric( vffirst[,locini:ncol( vffirst )] ) ) )
    colorm   <- vfcolormap( as.numeric( vfmiddlep[,locini:ncol( vfmiddle )] ) )
    vfmiddle <- as.character( round( as.numeric( vfmiddle[,locini:ncol( vfmiddle )] ) ) )
    colorl   <- vfcolormap( as.numeric( vflastp[,locini:ncol( vflast )] ) )
    vflast   <- as.character( round( as.numeric( vflast[,locini:ncol( vflast )] ) ) )
  }

  # add void points if not in the test of locations
  for( i in 1:nrow( bsxy ) ) {
    if( xmin < bsxy$x[i] & xmax > bsxy$x[i] & ymin < bsxy$y[i] & ymax > bsxy$y[i] ) {
      idx <- which( patternMap$xod == bsxy$x[i] & patternMap$yod == bsxy$y[i] )
      if( length( idx ) > 0 ) {
        colorf[idx,] <- c( 0, 0, 0 )
        colorm[idx,] <- c( 0, 0, 0 )
        colorl[idx,] <- c( 0, 0, 0 )
      } else {
        patternMap <- rbind( patternMap, c( bsxy$x[i], bsxy$y[i] ) )
        vffirst    <- c( vffirst,  NA )
        vfmiddle   <- c( vfmiddle,  NA )
        vflast     <- c( vflast,  NA )
        colorf     <- rbind( colorf, c( 0, 0, 0 ) )
        colorm     <- rbind( colorm, c( 0, 0, 0 ) )
        colorl     <- rbind( colorl, c( 0, 0, 0 ) )
      }
    }
  }

  if( vf$seye[1] == "OS" ) {
    xmin2 <- xmin
    xmin  <- -xmax
    xmax  <- -xmin2
    patternMap$xod <- -patternMap$xod
  }
  # get the Voronoi tesselation tiles tiles
  vftess  <- vftessellation( patternMap, dist = 3 )
  vftiles <- tile.list( vftess[[1]] )
  vfhull  <- vftess[[2]]

  ######################################################
  # first plot all graphs
  ######################################################
  # plot linear regression
  if( plotType == "vf" ) {
    ylab  <- "sensitivities"
    index <- vfindices$msens
  }
  if( plotType == "td" ) {
    ylab <- "MD"
    index <- vfindices$mtdev
  }
  if( plotType == "pd" ) {
    ylab <- "MPD 85th perc"
    index <- vfindices$mpdev
  }
  if( plotType == "pdghr" ) {
    ylab <- "MPD ranked"
    index <- vfindices$mpdev
  }
  opar <- par( no.readonly = TRUE )
  # PoPLR pointwise plot
  par( fig = c( 0.05, 0.45, 0.30, 0.65 ) )
  vfplot_plr( pres$sl, pres$pvall, pres$vfdata, colorMapType = colorMapType, colorScale = colorScale,
              showaxis = showaxis, colaxis = colaxis )
  # color scale
  par( new = TRUE )
  par( fig = c( 0.1, 0.325, 0.30, 0.33 ) )
  colormapgraph( ncol = 5, mapval = colorScale, notSeenAsBlack = FALSE )
  # OLS progression
  par( new = TRUE )
  par( fig = c( 0.00, 0.40, 0.02, 0.30 ) )
  progols( vfindices$tdate, index, ylab = ylab, txtfont = txtfont, pointsize = pointsize )
  # sensitivity plot first n visits
  par( new = TRUE )
  par( fig = c( 0.525, 1.00, 0.60, 0.95 ) )
  vfplotloc( vffirst, patternMap, vftiles = vftiles, vfhull = vfhull, loccol = colorf,
             xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
             showaxis = showaxis, colaxis = colaxis )
  # sensitivity plot middle n visits
  par( new = TRUE )
  par( fig = c( 0.525, 1.00, 0.30, 0.65 ) )
  vfplotloc( vfmiddle, patternMap, vftiles = vftiles, vfhull = vfhull, loccol = colorm,
             xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
             showaxis = showaxis, colaxis = colaxis )
  # sensitivity plot middle n visits
  par( new = TRUE )
  par( fig = c( 0.525, 1.00, 0.00, 0.35 ) )
  vfplotloc( vflast, patternMap, vftiles = vftiles, vfhull = vfhull, loccol = colorl,
             xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
             showaxis = showaxis, colaxis = colaxis )
  par( opar )

  ######################################################
  # create the text elements in the printouts
  ######################################################
  mainInfo   <- createviewport( "mainInfo",   left =  0.00, top =  0.00, width = 4.75, height = 0.40, pheight = mheight, pwidth = mwidth )
  infobox2   <- createviewport( "infobox2",   left =  6.40, top =  0.00, width = 1.40, height = 0.40, pheight = mheight, pwidth = mwidth )
  infobox3   <- createviewport( "infobox3",   left =  3.40, top = 10.90, width = 1.40, height = 0.30, pheight = mheight, pwidth = mwidth )
  textpoplr  <- createviewport( "textpoplr",  left =  0.20, top =  4.05, width = 3.40, height = 0.30, pheight = mheight, pwidth = mwidth )
  textvisitf <- createviewport( "textvisitf", left =  4.20, top =  0.55, width = 3.40, height = 0.30, pheight = mheight, pwidth = mwidth )
  textvisitm <- createviewport( "textvisitm", left =  4.20, top =  4.05, width = 3.40, height = 0.30, pheight = mheight, pwidth = mwidth )
  textvisitl <- createviewport( "textvisitl", left =  4.20, top =  7.55, width = 3.40, height = 0.30, pheight = mheight, pwidth = mwidth )
  results1   <- createviewport( "results1",   left =  0.25, top =  1.00, width = 0.50, height = 1.10, pheight = mheight, pwidth = mwidth )
  results2   <- createviewport( "results2",   left =  1.80, top =  1.00, width = 0.70, height = 1.10, pheight = mheight, pwidth = mwidth )
  results3   <- createviewport( "results3",   left =  2.60, top =  1.00, width = 1.00, height = 1.10, pheight = mheight, pwidth = mwidth )

  # create the list and then generate the tree and "push" it
  list <- vpList( mainInfo, infobox2, infobox3, textpoplr, textvisitf, textvisitm, textvisitl, results1, results2, results3 )
  tree <- vpTree( printout, list )
  
  pushViewport( tree )
  
  ######################################################
  # perimetry information
  ######################################################
  seekViewport( "mainInfo" )
  text <- vf$tperimetry[1]
  if( text == "sap" ) {
    text = "Static Automated Perimetry."
  } else if( text == "fdp" ) {
    text = "Frequency-doubling Perimetry."
  } else if( text == "csp" ) {
    text = "Contrast sensitivity Perimetry."
  }
  text <- paste( text, "Progression analysis", sep = " " )
  # ID
  text <- paste( text, "Subject ID: ", sep = "\n" )
  text <- paste( text, vf$id[1], ",", sep = "" )
  # age
  text <- paste( text, " age: from ", round( vf$sage[1] ), " to ", round( vf$sage[nrow( vf )] ), ",", sep = "" )
  # eye
  texteye <- paste( "eye:", vf$seye[1], sep = " " )
  if( vf$seye[1] == "OD" ) {
    texteye <- paste( texteye, "(right)", sep = " " )
  } else if ( vf$seye[1] == "OS" ) {
    texteye <- paste( texteye, "(left)", sep = " " )
  } else {
    texteye <- paste( texteye, "(which?)", sep = " " )
  }
  text <- paste( text, texteye, sep = " " )
  grid.text( text, x = 0.0, y = 1.0, just = c( "left", "top" ), gp = gpar( fontfamily = txtfont, fontsize = pointsize, fontface = "bold" ) )
  
  ######################################################
  # Details about printouts
  ######################################################
  seekViewport( "infobox2" )
  
  if( vf$tpattern[1] == "p24d2" ) {
    textpattern <- "Central 24-2"
  } else if( vf$tpattern[1] == "p30d2" ) {
    textpattern <- "Central 30-2"
  } else if( vf$tpattern[1] == "p10d2" ) {
    textpattern <- "Central 10-2"
  } else if( vf$tpattern[1] == "rnfl" ) {
    textpattern <- "RNFL-57"
  } else {
    textpattern <- "Unknown"
  }
  # algorithm
  if( vf$talgorithm[1] == "sitas" ) {
    textalgorithm <- "SITA standard"
  } else if( vf$talgorithm[1] == "sitaf" ) {
    textalgorithm <- "SITA fast"
  } else if( vf$talgorithm[1] == "fullt" ) {
    textalgorithm <- "Full threshold"
  } else {
    textalgorithm <- "Unknown"
  }
  
  text <- paste( textpattern, textalgorithm, sep = "\n" )
  grid.text( text, x = 1.00, y = 1.00, just = c( "right", "top" ), gp = gpar( fontfamily = txtfont, fontsize = pointsize ) )

  seekViewport( "textpoplr" )
  text <- "pointwise linear regression"
  if( plotType == "vf" )    text <- paste( text, "sensitivities" )
  if( plotType == "td" )    text <- paste( text, "TD" )
  if( plotType == "pd" )    text <- paste( text, "PD" )
  if( plotType == "pdghr" ) text <- paste( text, "PD ghr" )
  grid.text( text, x = 0.50, y = 0.50, just = c( "center", "center" ), gp = gpar( fontfamily = txtfont, fontsize = pointsize ) )

  ######################################################
  # Details about first, middle and last visits
  ######################################################
  seekViewport( "textvisitf" )

  if( grp == 1 ) {
    text <- paste( "first exam. Date: ", substr( tdatefirst, 6, 7 ), "/", substr( tdatefirst, 9, 10 ), "/", substr( tdatefirst, 1, 4 ), sep = "" )
  } else {
    text <- paste( "first ", grp, " exams. From ", substr( tdatefirst[1], 6, 7 ), "/", substr( tdatefirst[1], 9, 10 ), "/", substr( tdatefirst[1], 1, 4 ), " to ", substr( tdatefirst[grp], 6, 7 ), "/", substr( tdatefirst[grp], 9, 10 ), "/", substr( tdatefirst[grp], 1, 4 ), sep = "" )
  }
  grid.text( text, x = 0.50, y = 0.50, just = c( "center", "center" ), gp = gpar( fontfamily = txtfont, fontsize = pointsize ) )
  
  seekViewport( "textvisitm" )
  
  if( grp == 1 ) {
    text <- paste( "middle exam. Date: ", substr( tdatemiddle, 6, 7 ), "/", substr( tdatemiddle, 9, 10 ), "/", substr( tdatemiddle, 1, 4 ), sep = "" )
  } else {
    text <- paste( "middle ", grp, " exams. From ", substr( tdatemiddle[1], 6, 7 ), "/", substr( tdatemiddle[1], 9, 10 ), "/", substr( tdatemiddle[1], 1, 4 ), " to ", substr( tdatemiddle[grp], 6, 7 ), "/", substr( tdatemiddle[grp], 9, 10 ), "/", substr( tdatemiddle[grp], 1, 4 ), sep = "" )
  }
  grid.text( text, x = 0.50, y = 0.50, just = c( "center", "center" ), gp = gpar( fontfamily = txtfont, fontsize = pointsize ) )

  seekViewport( "textvisitl" )
  
  if( grp == 1 ) {
    text <- paste( "last exam. Date: ", substr( tdatelast, 6, 7 ), "/", substr( tdatelast, 9, 10 ), "/", substr( tdatelast, 1, 4 ), sep = "" )
  } else {
    text <- paste( "last ", grp, " exams. From ", substr( tdatelast[1], 6, 7 ), "/", substr( tdatelast[1], 9, 10 ), "/", substr( tdatelast[1], 1, 4 ), " to ", substr( tdatelast[grp], 6, 7 ), "/", substr( tdatelast[grp], 9, 10 ), "/", substr( tdatelast[grp], 1, 4 ), sep = "" )
  }
  grid.text( text, x = 0.50, y = 0.50, just = c( "center", "center" ), gp = gpar( fontfamily = txtfont, fontsize = pointsize ) )

  ######################################################
  # Overall statistics and results
  ######################################################
  lt      <- as.POSIXlt( vfindices$tdate )
  mon     <- lt$year * 12 + lt$mon
  yeardif <- ( mon - mon[1] ) / 12
  slm     <- summary( lm( index ~ yeardif ) )

  if( plotType == "vf" )    textindex <- "MS"
  if( plotType == "td" )    textindex <- "MD"
  if( plotType == "pd" )    textindex <- "PMD"
  if( plotType == "pdghr" ) textindex <- "PMDr"
  
  seekViewport( "results1" )
  
  text <- "Slope (per decade):"
  text <- paste( text, "PoPLR analysis p-value:", sep = "\n" )
  text <- paste( text, paste( textindex, "first",  grp, "visits:" ), sep = "\n" )
  text <- paste( text, paste( textindex, "middle", grp, "visits:" ), sep = "\n" )
  text <- paste( text, paste( textindex, "last",   grp, "visits:" ), sep = "\n" )
  text <- paste( text, paste( textindex, "projection (5 years):" ), sep = "\n" )
  grid.text( text, x = 0.00, y = 1.00, just = c( "left", "top" ), gp = gpar( fontfamily = txtfont, fontsize = pointsize ) )
  
  seekViewport( "results2" )
  
  text <- round( 10 * slm$coefficients[2,1], 1 )
  text <- paste( text, round( pres$pvall, 3 ), sep = "\n" )
  text <- paste( text, round( mean( index[idxfirst]  ), 2 ), sep = "\n" )
  text <- paste( text, round( mean( index[idxmiddle] ), 2 ), sep = "\n" )
  text <- paste( text, round( mean( index[idxlast]   ), 2 ), sep = "\n" )
  text <- paste( text, round( slm$coefficients[1,1] + slm$coefficients[2,1] * ( yeardif[length( yeardif )] + 5 ), 2 ), sep = "\n" )
  grid.text( text, x = 1.00, y = 1.00, just = c( "right", "top" ), gp = gpar( fontfamily = txtfont, fontsize = pointsize ) )
  
  seekViewport( "results3" )
  
  text  <- paste( "(p = ", round( slm$coefficients[2,4] / 2, 3 ), ")", sep = "" )
  grid.text( text, x = 0.00, y = 1.00, just = c( "left", "top" ), gp = gpar( fontfamily = txtfont, fontsize = pointsize ) )

  ######################################################
  # Details about printouts
  ######################################################
  seekViewport( "infobox3" )
  
  text <- paste( "norm vals: ", nv$nvname, sep = "" )
  text <- paste( text, ", visualFields ", packageDescription( "visualFields" )$Version, sep = "" )
  text <- paste( text, " (", substr( packageDescription( "visualFields" )$Date, 1, 4 ), ")", sep = "" )
  grid.text( text, x = 0.50, y = 0.00, just = c( "center", "bottom" ), gp = gpar( fontfamily = txtfont, fontsize = pointsize ) )

  # only if in save mode, then set device to off
  if( !is.null( filename ) ) {
    dev.off()
  }
}